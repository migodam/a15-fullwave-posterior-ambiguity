"""Full-wave multiple-cylinder E_z expansion, mu=1, outgoing Hankel-1.
Two finite circles in an infinite background; no Born/coupled-mode approximation.
Current coordinates are orthonormal in the volume L2 metric. Derived from exact
single-cylinder continuity and Graf translation. Cross-check against reused VIE.
"""
from functools import lru_cache
import numpy as np
from scipy.special import jv,jvp,hankel1,h1vp
from scipy.linalg import solve
from numpy.polynomial.legendre import leggauss
from physics import Family,PhysicalModel

class CylinderModel:
    def __init__(self,family=None,order=12,radial_order=48):
        self.family=family or Family();self.order=order;self.m=np.arange(-order,order+1);self.d=len(self.m)
        self.receivers=PhysicalModel(16,self.family).receivers
        x,w=leggauss(radial_order);self.radial=self.family.radius*(x+1)/2;self.rw=w*self.family.radius/2*self.radial
        self.solve_count=0

    def components(self,theta,k,coupling=True):
        s,delta,gamma=np.asarray(theta,float);f=self.family;r=f.radius
        if not (f.gap_min-1e-12<=s<=f.gap_max+1e-12 and gamma>=0 and f.epsilon_base+delta>1):raise ValueError('invalid physical theta')
        centers=np.array([[f.target_x-2*r-s,0],[f.target_x,0]])
        eps=np.array([f.epsilon_base+1j*f.bright_gamma*f.k_ref/k,f.epsilon_base+delta+1j*gamma*f.k_ref/k])
        ki=k*np.sqrt(eps);m=self.m;x=k*r
        parts=[self.body_coefficients(float(k),complex(e)) for e in eps]
        T,D,integ,norm=map(np.array,zip(*parts))
        inject=norm*D
        radiation=T/inject
        U=self.translation(float(k),float(s)) if coupling else np.zeros((2*self.d,2*self.d),complex)
        M=np.eye(2*self.d)-inject.ravel()[:,None]*U*radiation.ravel()[None,:]
        incident=np.exp(1j*k*centers[:,0,None])*(1j**m)[None,:]
        rhs=(inject*incident).ravel()
        return dict(M=M,rhs=rhs,centers=centers,epsilon=eps,ki=ki,T=T,D=D,integrals=integ,norm=norm,inject=inject,radiation=radiation,U=U)

    @lru_cache(maxsize=8192)
    def body_coefficients(self,k,eps):
        m=self.m;r=self.family.radius;ki=k*np.sqrt(eps);x=k*r
        J=jv(m,x);Jp=jvp(m,x);H=hankel1(m,x);Hp=h1vp(m,x)
        Z=jv(m,ki*r);Zp=jvp(m,ki*r)
        T=(ki*Zp*J-k*Jp*Z)/(k*Hp*Z-ki*Zp*H)
        D=(J+T*H)/Z
        integ=np.sum(abs(jv(m[None,:],ki*self.radial[:,None]))**2*self.rw[:,None],axis=0)*2*np.pi
        norm=(eps-1)*np.sqrt(integ)
        return T,D,integ,norm

    @lru_cache(maxsize=512)
    def translation(self,k,gap):
        U=np.zeros((2*self.d,2*self.d),complex);mn=self.m[None,:]-self.m[:,None]
        for a,b,phi in [(0,1,np.pi),(1,0,0.)]:
            U[a*self.d:(a+1)*self.d,b*self.d:(b+1)*self.d]=hankel1(mn,k*(2*self.family.radius+gap))*np.exp(1j*mn*phi)
        return U

    @lru_cache(maxsize=512)
    def receiver_basis(self,k,gap):
        f=self.family;centers=[[f.target_x-2*f.radius-gap,0],[f.target_x,0]];blocks=[]
        for center in centers:
            xy=self.receivers-np.asarray(center);rr=np.linalg.norm(xy,axis=1);angle=np.arctan2(xy[:,1],xy[:,0])
            blocks.append(hankel1(self.m[None,:],k*rr[:,None])*np.exp(1j*angle[:,None]*self.m))
        return np.concatenate(blocks,axis=1)

    def external_operator(self,points,k,state):
        if points is self.receivers:
            gap=float(state['centers'][1,0]-state['centers'][0,0]-2*self.family.radius)
            # Round only caching key to remove center subtraction roundoff, not material parameters.
            return self.receiver_basis(float(k),round(gap,14))*state['radiation'].ravel()[None,:]
        pts=np.asarray(points);blocks=[]
        for b in range(2):
            xy=pts-state['centers'][b];rr=np.linalg.norm(xy,axis=1);angle=np.arctan2(xy[:,1],xy[:,0])
            blocks.append(hankel1(self.m[None,:],k*rr[:,None])*np.exp(1j*angle[:,None]*self.m)*state['radiation'][b])
        return np.concatenate(blocks,axis=1)

    def solve(self,theta,k,acquisition=None,coupling=True,fields=False):
        st=self.components(theta,k,coupling);z=solve(st['M'],st['rhs'],check_finite=False)
        residual=float(np.linalg.norm(st['M']@z-st['rhs'])/max(np.linalg.norm(st['rhs']),1e-30))
        energy=np.sum(abs(z.reshape(2,self.d)/(st['epsilon'][:,None]-1))**2,axis=1)
        q=float(np.log(energy[1]/(np.pi*self.family.radius**2)))
        S=self.external_operator(self.receivers,k,st);y=S@z
        if not np.isfinite(q) or not np.isfinite(y).all() or residual>1e-10:raise FloatingPointError(f'residual={residual}')
        out={'theta':list(map(float,theta)),'k':float(k),'Q':q,'y':y if acquisition is None else y[acquisition],
             'absorption':float(k*np.sum(st['epsilon'].imag*energy)),
             'current_energy':float(np.sum(abs(z[self.d:])**2)),
             'energy_by_body':energy,'residual':residual,'active_cells':2*self.d,'coupling':bool(coupling)}
        self.solve_count+=1
        if fields:out.update(state=st,z=z,S=S,internal_coefficients=z.reshape(2,self.d)/st['norm'])
        return out

    def field(self,solution,points):
        st=solution['state'];k=solution['k'];pts=np.asarray(points)
        out=np.exp(1j*k*pts[:,0]);outside=np.ones(len(pts),bool)
        for b in range(2):
            xy=pts-st['centers'][b];rr=np.linalg.norm(xy,axis=1);phi=np.arctan2(xy[:,1],xy[:,0]);inside=rr<self.family.radius
            outside &= ~inside
            out[inside]=np.sum(solution['internal_coefficients'][b]*jv(self.m[None,:],st['ki'][b]*rr[inside,None])*np.exp(1j*phi[inside,None]*self.m),axis=1)
        if outside.any():out[outside]+=self.external_operator(pts[outside],k,st)@solution['z']
        return out

    def flux_diagnostics(self,theta,k,n_angles=256):
        sol=self.solve(theta,k,fields=True);a=np.arange(n_angles)*2*np.pi/n_angles;normal=np.c_[np.cos(a),np.sin(a)];pts=2*normal
        st=sol['state'];es=np.zeros(n_angles,complex);ds=es.copy()
        for b in range(2):
            xy=pts-st['centers'][b];rr=np.linalg.norm(xy,axis=1);phi=np.arctan2(xy[:,1],xy[:,0]);rad=xy/rr[:,None];tang=np.c_[-rad[:,1],rad[:,0]]
            phase=np.exp(1j*phi[:,None]*self.m);bb=st['radiation'][b]*sol['z'][b*self.d:(b+1)*self.d]
            H=hankel1(self.m[None,:],k*rr[:,None]);Hp=h1vp(self.m[None,:],k*rr[:,None]);es+=np.sum(H*phase*bb,axis=1)
            ds+=np.sum(phase*bb*(k*Hp*np.sum(rad*normal,axis=1)[:,None]+H*1j*self.m/rr[:,None]*np.sum(tang*normal,axis=1)[:,None]),axis=1)
        ei=np.exp(1j*k*pts[:,0]);di=1j*k*normal[:,0]*ei;measure=4*np.pi/n_angles
        sca=float(measure/k*np.imag(np.vdot(es,ds)));ext=float(-measure/k*np.imag(np.vdot(ei,ds)+np.vdot(es,di)));bal=ext-sca-sol['absorption']
        return {'scattered':sca,'extinction':ext,'absorption':sol['absorption'],'balance_relative':abs(bal)/max(abs(ext),1e-12)}

    def parameter_jacobian(self,theta,k,step=1e-5):
        """Differentiate assembled full-wave equation, including moving centers.
        The mode metric depends on theta: both matrix and RHS must differentiate.
        Q includes the derivative of the current-to-field normalization.
        """
        th=np.array(theta,float);sol=self.solve(th,k,fields=True);st=sol['state'];z=sol['z'];dy=[];dq=[];dz=[]
        for i,h in enumerate(step*np.array([.1,10.,.1])):
            tp=th.copy();tm=th.copy();tp[i]+=h;tm[i]-=h
            sp=self.components(tp,k);sm=self.components(tm,k)
            dr=(sp['rhs']-sm['rhs'])/(2*h)-(sp['M']-sm['M'])@z/(2*h)
            zz=solve(st['M'],dr,check_finite=False);dz.append(zz)
            ss=(self.external_operator(self.receivers,k,sp)-self.external_operator(self.receivers,k,sm))/(2*h)
            dy.append(sol['S']@zz+ss@z)
            norm=st['epsilon'][1]-1;dn=(sp['epsilon'][1]-sm['epsilon'][1])/(2*h)
            e=z[self.d:]/norm;de=zz[self.d:]/norm-z[self.d:]*dn/norm**2
            dq.append(2*np.real(np.vdot(e,de))/np.vdot(e,e).real)
        return {'y':np.array(dy).T,'Q':np.array(dq),'z':np.array(dz).T}

    def parameter_jvp(self,theta,k,v,task='y'):return self.parameter_jacobian(theta,k)[task]@np.asarray(v)
    def parameter_vjp(self,theta,k,cotangent,task='y'):return np.real(self.parameter_jacobian(theta,k)[task].conj().T@np.asarray(cotangent))
