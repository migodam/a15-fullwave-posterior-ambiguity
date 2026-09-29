"""Passive scalar E_z volume-integral model. e^{-i omega t}, epsilon'' >= 0.
J=chi E is contrast source, NOT physical heating. Spatial geometry is continuous.
The A14 kernel and corrected equal-area self-cell are reused verbatim.
"""
from dataclasses import dataclass
from functools import lru_cache
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv, jvp as bessel_jvp, h1vp
from physics_backend import make_grid, green_matrix, self_cell_green


def pack(z):
    z=np.asarray(z).ravel();return np.r_[z.real,z.imag]


def circle_fractions(points,h,center,radius):
    """Exact disk-square overlap via segmented analytic sqrt antiderivative."""
    x=points[:,0]-center[0];y=points[:,1]-center[1];r=radius
    nearest=np.maximum(np.abs(x)-h/2,0)**2+np.maximum(np.abs(y)-h/2,0)**2
    farthest=(np.abs(x)+h/2)**2+(np.abs(y)+h/2)**2
    out=np.zeros(len(x));out[farthest<=r*r]=1
    boundary=np.flatnonzero((nearest<r*r)&(farthest>r*r))
    def prim(t):
        t=np.clip(t,-r,r)
        return .5*(t*np.sqrt(max(0.,r*r-t*t))+r*r*np.arcsin(t/r))
    for i in boundary:
        xa=max(-r,x[i]-h/2);xb=min(r,x[i]+h/2);ya=y[i]-h/2;yb=y[i]+h/2
        breaks=[xa,xb]
        for yy in [ya,yb]:
            if abs(yy)<r:
                t=np.sqrt(r*r-yy*yy)
                breaks += [v for v in [-t,t] if xa<v<xb]
        breaks=sorted(set(breaks));area=0.
        for a,b in zip(breaks[:-1],breaks[1:]):
            z=np.sqrt(max(0.,r*r-((a+b)/2)**2))
            if min(yb,z)<=max(ya,-z):continue
            c=(yb if yb<z else 0.)-(ya if ya>-z else 0.)
            m=(0 if yb<z else 1)+(0 if ya>-z else 1)
            area+=c*(b-a)+m*(prim(b)-prim(a))
        out[i]=np.clip(area/h**2,0,1)
    return out


@dataclass(frozen=True)
class Family:
    radius:float=.12
    epsilon_base:float=12.
    target_x:float=.18
    bright_gamma:float=.002
    k_ref:float=2*np.pi
    gap_min:float=.02
    gap_max:float=.16


class PhysicalModel:
    def __init__(self,n=48,family=None,backend='cpu',threads=1):
        self.n=n;self.family=family or Family();self.backend=backend
        self.all_points,self.h=make_grid(n);f=self.family;r=f.radius
        # Fixed union covers every translated bright disk and the fixed target.
        x=self.all_points[:,0];y=self.all_points[:,1]
        lo=f.target_x-2*r-f.gap_max;hi=f.target_x-2*r-f.gap_min
        dx=np.maximum(np.maximum(lo-x,x-hi),0)
        bright=(np.maximum(dx-self.h/2,0)**2+np.maximum(abs(y)-self.h/2,0)**2)<r*r
        dark=circle_fractions(self.all_points,self.h,(f.target_x,0),r)>0
        self.indices=np.flatnonzero(bright|dark);self.points=self.all_points[self.indices]
        self.roi=circle_fractions(self.points,self.h,(f.target_x,0),r)
        # Unique non-mirror external channels in upper semicircle; near field last.
        angles=np.array([0,np.pi,np.pi/4,3*np.pi/4,np.pi/8,3*np.pi/8,5*np.pi/8,7*np.pi/8,
                         np.pi/16,3*np.pi/16,5*np.pi/16,7*np.pi/16,9*np.pi/16,11*np.pi/16,13*np.pi/16,15*np.pi/16])
        self.receivers=np.vstack([2*np.c_[np.cos(angles),np.sin(angles)],
                                [f.target_x+r+.04,0],[f.target_x,r+.04]])
        self._operators={};self.solve_count=0

    def fractions(self,gap):
        f=self.family
        return circle_fractions(self.points,self.h,(f.target_x-2*f.radius-gap,0),f.radius),self.roi

    def material(self,theta,k):
        s,d,g=np.asarray(theta,float);f=self.family
        if not (f.gap_min-1e-12<=s<=f.gap_max+1e-12 and g>=0 and f.epsilon_base+d>1 and k>0):
            raise ValueError('nonphysical or out-of-geometry-support parameters')
        b,t=self.fractions(s)
        chi=b*(f.epsilon_base-1+1j*f.bright_gamma*f.k_ref/k)+t*(f.epsilon_base+d-1+1j*g*f.k_ref/k)
        return {'chi':chi,'bright_fraction':b,'target_fraction':t,'roi':self.roi,'theta':np.asarray(theta,float)}

    def operators(self,k):
        key=float(k)
        if key not in self._operators:
            gd=k*k*self.h*self.h*green_matrix(self.points,self.points,k)
            np.fill_diagonal(gd,k*k*self_cell_green(k,self.h))
            gs=k*k*self.h*self.h*green_matrix(self.points,self.receivers,k)
            inc=np.exp(1j*k*self.points[:,0])
            if len(self._operators)>=3:self._operators.pop(next(iter(self._operators)))
            self._operators[key]=(gd,gs,inc)
        return self._operators[key]

    def solve(self,theta,k,acquisition=None,coupling=True,fields=False):
        started=time.perf_counter();m=self.material(theta,k);chi=m['chi'];gd,gs,inc=self.operators(k)
        active=np.flatnonzero(abs(chi)>0)
        g=gd[np.ix_(active,active)].copy()
        labels=m['target_fraction'][active]>0
        if not coupling:g[labels[:,None]!=labels[None,:]]=0
        mat=np.eye(len(active),dtype=complex)-chi[active,None]*g;rhs=chi[active]*inc[active]
        if self.backend=='cuda':
            import torch
            a=torch.as_tensor(mat,device='cuda');b=torch.as_tensor(rhs,device='cuda')
            ja=torch.linalg.solve(a,b).cpu().numpy();lu=None
        else:
            lu=lu_factor(mat,check_finite=False);ja=lu_solve(lu,rhs,check_finite=False)
        j=np.zeros(len(chi),complex);j[active]=ja
        e=inc+gd@j
        if not coupling:e[active]=inc[active]+g@ja
        y=gs@j
        residual=float(np.linalg.norm(mat@ja-rhs)/max(np.linalg.norm(rhs),1e-30))
        if not np.isfinite(e).all() or residual>1e-10:raise FloatingPointError(f'invalid solve residual {residual}')
        q=float(np.log(np.sum(self.roi*abs(e)**2)/self.roi.sum()))
        # Cross section normalized by incident intensity |E0|^2/(2 Z0), c=Z0=1.
        absorption=float(k*self.h**2*np.sum(chi.imag*abs(e)**2))
        self.solve_count+=1
        result={'theta':list(map(float,theta)),'k':float(k),'Q':q,'y':y if acquisition is None else y[acquisition],
                'absorption':absorption,'current_energy':float(self.h**2*np.sum(self.roi*abs(j)**2)),
                'max_field':float(abs(e).max()),'residual':residual,'seconds':time.perf_counter()-started,
                'active_cells':len(active),'coupling':bool(coupling)}
        if fields:result.update(J=j,E=e,chi=chi,incident=inc,active=active,matrix=mat,lu=lu,GD=gd,GS=gs)
        return result

    def material_jacobian(self,theta,k,step=1e-6):
        theta=np.asarray(theta,float);f=self.family
        b,t=self.fractions(theta[0]);bp,_=self.fractions(theta[0]+step);bm,_=self.fractions(theta[0]-step)
        db=(bp-bm)/(2*step)
        return np.c_[db*(f.epsilon_base-1+1j*f.bright_gamma*f.k_ref/k),t,1j*f.k_ref/k*t]

    def parameter_jacobian(self,theta,k):
        f=self.solve(theta,k,fields=True);dc=self.material_jacobian(theta,k)
        # Union-space tangent includes any zero-chi cells with nonzero geometry derivative.
        mat=np.eye(len(f['chi']),dtype=complex)-f['chi'][:,None]*f['GD']
        dj=lu_solve(lu_factor(mat),f['E'][:,None]*dc)
        de=f['GD']@dj;dy=f['GS']@dj
        dq=2*np.real((self.roi*np.conj(f['E']))@de)/np.sum(self.roi*abs(f['E'])**2)
        return {'J':dj,'E':de,'y':dy,'Q':dq}

    def parameter_jvp(self,theta,k,direction,task='y'):
        return self.parameter_jacobian(theta,k)[task]@np.asarray(direction,float)

    def parameter_vjp(self,theta,k,cotangent,task='y'):
        jac=self.parameter_jacobian(theta,k)[task]
        return np.real(np.conj(jac).T@np.asarray(cotangent))

    def flux_diagnostics(self,theta,k,n_angles=512):
        f=self.solve(theta,k,fields=True);angles=np.arange(n_angles)*2*np.pi/n_angles
        pts=2*np.c_[np.cos(angles),np.sin(angles)];normals=pts/2
        diff=pts[:,None,:]-self.points[None,:,:];dist=np.linalg.norm(diff,axis=-1)
        gs=k*k*self.h**2*(1j/4)*hankel1(0,k*dist)
        gn=k*k*self.h**2*(-1j*k/4)*hankel1(1,k*dist)*np.einsum('rpi,ri->rp',diff,normals)/dist
        es=gs@f['J'];ds=gn@f['J'];ei=np.exp(1j*k*pts[:,0]);di=1j*k*normals[:,0]*ei
        measure=4*np.pi/n_angles
        scattered=float(measure/k*np.imag(np.vdot(es,ds)))
        extinction=float(-measure/k*np.imag(np.vdot(ei,ds)+np.vdot(es,di)))
        balance=extinction-scattered-f['absorption']
        return {'scattered_cross_section':scattered,'extinction_cross_section':extinction,
                'absorption_cross_section':f['absorption'],'balance_abs':balance,
                'balance_relative':abs(balance)/max(abs(extinction),abs(scattered)+abs(f['absorption']),1e-12)}


def cylinder_exact(points,k,radius,epsilon,center=(0.,0.),order=32):
    """Scalar E_z cylinder: continuity of E and its normal derivative (mu=1)."""
    points=np.asarray(points);xy=points-np.asarray(center);rr=np.linalg.norm(xy,axis=1);ang=np.arctan2(xy[:,1],xy[:,0])
    m=np.arange(-order,order+1);ki=k*np.sqrt(complex(epsilon));x=k*radius;z=ki*radius
    jx=jv(m,x);jx1=bessel_jvp(m,x);jz=jv(m,z);jz1=bessel_jvp(m,z)
    hx=hankel1(m,x);hx1=h1vp(m,x)
    a=(ki*jz1*jx-k*jx1*jz)/(k*hx1*jz-ki*jz1*hx)
    coeff=1j**m*np.exp(1j*k*center[0]);out=np.zeros(len(points),complex)
    ext=rr>=radius
    out[ext]=np.sum(coeff[None,:]*a[None,:]*hankel1(m[None,:],k*rr[ext,None])*np.exp(1j*ang[ext,None]*m),axis=1)
    if (~ext).any():
        c=(jx+a*hx)/jz
        out[~ext]=np.sum(coeff[None,:]*c[None,:]*jv(m[None,:],ki*rr[~ext,None])*np.exp(1j*ang[~ext,None]*m),axis=1)
    return out


def physical_readout(theta,conditions,model):
    return model.solve(theta,conditions['k'],conditions.get('acquisition'),conditions.get('coupling',True),conditions.get('fields',False))
