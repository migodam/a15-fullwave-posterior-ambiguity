# Supplement to “Current-Space Enrichment in Electromagnetic Inversion: Material-Update Closure and Design Limits”

This supplement gives coordinate conventions, detailed local proofs, and the existing A12 execution evidence. The main manuscript has two principal theorems: paired step-change closure and unchanged-full-goal domination. Certificate and constraint results are supporting conditions. No experiment, retraining, reconstruction, or remote computation was added for this final writing. Reference numbers refer to the final manuscript.

## S1. Physical coordinates and the implemented vector-Maxwell model

All material normal equations are over real Hilbert spaces, with material coordinates \(s\in\mathbb R^d\). For complex material contrast, the real and imaginary grid values are stacked, and the real material expansion satisfies \(Q_T^\top M_\chi Q_T=I\). A complex-output map \(F\) of real material coordinates has real adjoint \(v\mapsto\operatorname{Re}(F^Hv)\); equivalently, use its fully realified matrix. For a complex current operator,

$$
\mathbb R(F)=\begin{bmatrix}\operatorname{Re}F&-\operatorname{Im}F\\
\operatorname{Im}F&\operatorname{Re}F\end{bmatrix}.
$$

For \(P\) illuminations, \(L\), \(S\), \(V\), and \(C\) lift blockwise, while the \(B_p\) stack vertically and share the same \(s\). In particular, \(C_{\rm blk}=I_P\otimes\mathbb R(C)\) has \(2Pr_c\) columns. One cannot choose separate material updates for separate illuminations, or replace the vertically stacked Jacobian covariance by a sum that omits cross-illumination blocks.

The current and material metrics are normalized before all orthogonalizations. Orthogonal rotations of an admissible physical tangent preserve the assertions. Different Gaussian and voxel spans are different admissible tangents. Their numerical comparison does not establish an invariant ranking across different physical spaces.

The actual reused DDA models an isotropic nonmagnetic dielectric, \(\chi=\epsilon_r-1\), on cell centers in a box of side 1.5. A cell has volume \(v=(1.5/n)^3\). With \(\boldsymbol n=\boldsymbol r/\|\boldsymbol r\|\), the off-diagonal electric dyadic is

$$
G(\boldsymbol r)=\frac{e^{ikr}}{4\pi r}
\left[k^2(I-\boldsymbol n\boldsymbol n^\top)
+\left(\frac{ik}{r}-\frac1{r^2}\right)(I-3\boldsymbol n\boldsymbol n^\top)\right]. \tag{S1}
$$

Self blocks are excluded from \(G_{\rm off}\). The implemented cell polarizability and derivative are

$$
a(\chi)=\frac{3v\chi}{\chi+3-i\,3c_kv\chi},\qquad
c_k=\frac{k^3}{6\pi},\qquad
a'(\chi)=\frac{9v}{(\chi+3-i\,3c_kv\chi)^2}. \tag{S2}
$$

For dipole moments \(p_p\), the normalized current is \(c_p=p_p/\sqrt v\). The state, total field, tangent, and measurement equations are

$$
\begin{aligned}
L&=I-D_aG_{\rm off},&
Lc_p&=D_ae_p^{\rm inc}/\sqrt v,\\
E_p&=e_p^{\rm inc}+G_{\rm off}(\sqrt v\,c_p),&
L\delta c_p&=D_{a'E_p}Q_Ts/\sqrt v,\\
y_p&=\sqrt v\,G_{\rm rec}c_p.
\end{aligned} \tag{S3}
$$

The material-to-vector map in the tangent multiplies each cell's scalar material perturbation by its three field components. It is not a scalar diagonal on an unrelated vector space. Polarizability is not contrast; dropping \(a'(\chi)\) changes the Jacobian. At zero contrast \(a'=v\), so zero contrast does not remove the material tangent. Equation (S2) incorporates Clausius–Mossotti and radiation-reaction terms in one expression; the final denominator, not an intermediate rearrangement at \(\chi=-3\), specifies its singularities. These formulas document the reused code. Reference [11] supplies DDA background, not an audited derivation of this exact implementation.

There are six transverse plane-wave illuminations and 64 receivers on a sphere of radius 5, with two transverse readout functionals each. The receiver map retains the dyadic near-field terms. Measured scattered fields and tangent actions share the scalar normalization \(\|y_{\rm clean}\|\); no estimated noise-covariance whitening was used in the experiments. This clean norm remains known and fixed in the synthetic noise checks.

The material metric is \(vI\) for each complex component. Gaussian columns are \(\exp[-\|x_i-\xi_\nu\|^2/(2\times0.28^2)]\), with \(\xi_\nu\in\{-0.45,0,0.45\}^3\). An SVD of the volume-weighted columns produces the orthonormal physical tangent. Voxel expansion is \(\delta\chi=(s_{\rm re}+is_{\rm im})/\sqrt v\). Smooth Gaussian profiles define the internal truth, with object-dependent weights and a constant background; their fixed-box support is truncated at the outer boundary. No sharp internal junctions were tested.

## S2. Feedback factorization, bounds, and actual-ROM drift

Let \(A=I-L\), \(V^*V=C^*C=I\), \(V^*C=0\), and \(D_V=I-V^*AV\). All block formulas first hold in complex current coordinates and then lift to the real shared-material system. The enlarged Galerkin equation is

$$
\begin{bmatrix}D_V&-A_{VC}\\-A_{CV}&I-A_{CC}\end{bmatrix}
\begin{bmatrix}u\\a\end{bmatrix}
=\begin{bmatrix}B_V\\B_C\end{bmatrix}. \tag{S4}
$$

Solving the first row yields \(u=D_V^{-1}B_V+D_V^{-1}A_{VC}a\). Substitution into the second row gives

$$
\Gamma_Ca=N_C,\quad
\Gamma_C=I-A_{CC}-A_{CV}D_V^{-1}A_{VC},\quad
N_C=B_C+A_{CV}D_V^{-1}B_V. \tag{S5}
$$

Thus \(K_C=K_0+\bar C\Gamma_C^{-1}N_C\), where \(\bar C=C+VD_V^{-1}A_{VC}\). Since \(R_0C=0\), \(\bar C=(I-R_0L)C\). Direct multiplication gives \(V^*L\bar C=0\) and \(C^*L\bar C=\Gamma_C\). Observation gives \(\Delta J=Y_C\Gamma_C^{-1}N_C\). The injection-side adjoint is

$$
N_C^*=B^*(I-R_0^*L^*)C. \tag{S6}
$$

These formulas require invertibility of the retained and Schur blocks, not self-adjointness or normality of \(L\). Near-singular blocks invalidate an unqualified numerical closure claim even if a computed bank remains usable as a proposal space.

If \(\|A_{VV}\|\le\eta<1\), the Neumann-series bound gives \(\|D_V^{-1}\|\le(1-\eta)^{-1}\). It follows that

$$
\|Y_C-SC\|\le\frac{\|SV\|\|A_{VC}\|}{1-\eta},\qquad
\|N_C-B_C\|\le\frac{\|A_{CV}\|\|B_V\|}{1-\eta}. \tag{S7}
$$

When \(\sigma_{\min}(SC)>0\), divide the first bound applied to a coefficient vector by \(\|SCa\|\ge\sigma_{\min}(SC)\|a\|\). If \(SC\) has dark directions, no such all-direction relative bound follows; it can only be stated on a specified restricted coefficient subspace with a positive lower singular value. For \(F_C=A_{CC}+A_{CV}D_V^{-1}A_{VC}\), the sufficient condition \(\|F_C\|<1\) gives \(\|\Gamma_C^{-1}\|\le(1-\|F_C\|)^{-1}\). These are sufficient conditions, not numerical classifications established for every A12 scene.

The full-minus-actual-ROM diagnostic uses a different comparison. Write \(P=VV^*\), \(Q=I-P\), and let \(R_P\) invert \(PLP\) on \(\operatorname{Range}P\), extended by zero. The projected full tangent equation is

$$
PK=R_PPB+R_PPAQ\,QK.
$$

Subtract \(\widehat K=R_PP\widehat B\) and add \(QK\):

$$
K-\widehat K=QK+R_PPAQ\,QK+R_PP(B-\widehat B). \tag{S8}
$$

Only retained-block invertibility is needed for this identity. Applying \(S\) gives the three recorded observable contributions. Their vectors may interfere. The last term vanishes when the ROM interpolates the full primal field at the anchor and forms injection consistently, as in the enlargement comparison; it is intentionally observable in the unanchored Gs24 diagnostic.

For a genuinely changing ROM field, \(B_i s=Db[s]-DL[s]j_i\) and \(\delta B s=-DL[s](j_1-j_0)\). The enlarged inverse obeys

$$
R_C=R_0+\bar C\Gamma_C^{-1}C^*(I-LR_0).
$$

Hence \(K_1-K_0\) includes both the frozen-injection correction and \(R_C\delta B\). Its rank need not satisfy the small-bank bound stated for fixed \(B\). If the residual or regularizer also changes, its additional normal-equation terms must be included. The main theorem deliberately freezes all these quantities.

## S3. Effective support, minimality, inertia, and degeneracies

Let \(D=J_C-J_0\), with \(J_0,J_C,\Lambda\) fixed, and \(\Delta H=J_C^*J_C-J_0^*J_0\). Both \(H_0\) and \(H_C\) are positive definite. For common \(r,\ell\), subtraction gives

$$
H_0(s_C-s_0)=-\Delta Hs_C-D^*r. \tag{S9}
$$

The identity \(\Delta H=J_0^*D+D^*J_C\) shows that \(\operatorname{Range}[\Delta H,D^*]\) is contained in \(\operatorname{Range}[J_0^*D,D^*]\). Rearranging to \(J_0^*D=\Delta H-D^*J_C\) proves the reverse inclusion. Multiplying by invertible \(M_0\) yields the effective support in main (13).

For the minimality statement, take any independent pair \((x,r)\) and choose the common prior-linear term \(\ell=-H_Cx-J_C^*r\). Then \(s_C=x\), and (S9) ranges over \(-M_0\operatorname{Range}[\Delta H,D^*]\). Conversely every step difference belongs to that range. Thus the span of all differences over arbitrary \(r\) and shared \(\ell\) is exactly \(\mathcal E_{\rm pert}\). The quantifier on \(\ell\) is essential: a fixed prior gradient and a single residual may activate a smaller space. This is not a claim that the implemented raw bank is always minimal.

Let \(D=\widehat Y\widehat T\) have full-column-rank \(\widehat Y\in\mathbb R^{m\times q}\) and full-row-rank \(\widehat T\in\mathbb R^{q\times d}\). Then

$$
\operatorname{Range}(J_0^*D)=\operatorname{Range}(J_0^*\widehat Y),\quad
\operatorname{Range}(D^*)=\operatorname{Range}(\widehat T^*).
$$

Put \(U=J_0^*\widehat Y\), \(W=[U,\widehat T^*]\), and \(F=\widehat Y^*\widehat Y\). Expansion gives

$$
\Delta H=W\begin{bmatrix}0&I\\I&F\end{bmatrix}W^*. \tag{S10}
$$

The middle quadratic form is \(2x^*y+y^*Fy\). With \(\widetilde x=x+Fy/2\), it becomes \(2\widetilde x^*y\), which has \(q\) positive and \(q\) negative directions. If \(W\) has full column rank, a thin QR factorization reduces its nonzero inertia to a nonsingular congruence of that middle matrix; therefore \(\Delta H\) has rank \(2q\) and inertia \((q,q,d-2q)\). If the two signature ranges overlap or have dependencies, only the upper bound \(2q\) applies. This is established matrix algebra used to interpret the electromagnetic factors, not a new inertia theorem.

For the main two-coordinate example, choose current dimension two, \(L=I\), \(V=(1,0)^\top\), \(C=(0,1)^\top\), \(S=[1,1]\), and \(B=I\). Then \(R_0=\operatorname{diag}(1,0)\), \(J_0=[1,0]\), \(J_C=[1,1]\), \(Y_C=1\), and \(N_C=[0,1]\). With \(\Lambda=I,r=1,\ell=0\), the two material lift ranges are the two coordinate axes and the step difference has both components. This is the exact analytic illustration supplied by the final theory report, not a newly measured physical counterexample.

The following degeneracies follow directly:

| Condition | Consequence |
|---|---|
| \(N_C=0\) | No new material excitation and \(D=0\). |
| \(Y_C=0\) | No dressed observation and \(D=0\). |
| \(J_0^*\widehat Y=0\) | Injection alone contains the effective support; \(\Delta H=D^*D\succeq0\). |
| \(\operatorname{Range}U\subseteq\operatorname{Range}\widehat T^*\) | Injection alone suffices. |
| \(\operatorname{Range}\widehat T^*\subseteq\operatorname{Range}U\) | Observation alone suffices. |
| \(\Delta Hs_0+D^*r=0\) | The particular step is unchanged, even if \(D\ne0\). |

For the last row, the condition makes \(s_0\) satisfy the enlarged normal equation; uniqueness gives \(s_C=s_0\). Necessity follows by subtracting the equations. With \(P\) illuminations and \(r_c\) shared complex current columns, \(q\le2Pr_c\) and \(\dim\mathcal E_{\rm pert}\le\min(d,4Pr_c)\). The A12 measured rank 48 is the raw paired-bank rank after numerical revelation, not a separately measured minimal effective rank.

## S4. Full-goal optimization and defect alignment

For an independent material basis \(Z\), expanding the unchanged full quadratic gives

$$
\Phi(s_0+Za)=\Phi(s_0)+a^*Z^*d_0+\frac12a^*Z^*HZa.
$$

Minimization yields the main repair formula. Since \(s_C\in s_0+\mathcal E_0(C)\), the specified direct step is a feasible competitor. Completing the square gives the full quadratic and energy ordering. Numerically discarding nonzero signatures may lose that inclusion; the measured closure residual is then a separate required diagnostic.

For \(U_H=H^{1/2}Z\), its Euclidean projector is \(P=U_H(U_H^*U_H)^{-1}U_H^*\). Inserting this projector gives main (19). Thus a structurally complete envelope with \(H^{-1/2}d_0\perp H^{1/2}\mathcal E\) yields zero gain. By contrast, for \(d_0\ne0\) and \(M_0\succ0\), a line search along \(-M_0d_0\) yields positive quadratic gain because \(d_0^*M_0d_0>0\). This conditional algebraic observation explains why structural closure alone cannot order all material spaces; it is not an unreported experiment or a claim that every Krylov bank wins universally.

For redundant columns the inverse is replaced by a Moore–Penrose inverse on the effective range. If earlier paid directions \(U\) have already been optimized, \(U^*d=0\), remove overlap with \(D_Z=(I-\Pi_U^H)Z\), where \(\Pi_U^H=U(U^*HU)^{-1}U^*H\). Then

$$
\delta s=-D_ZG^\dagger b,\quad
G=D_Z^*HD_Z,\quad b=D_Z^*d,\quad
\Delta\Phi=\frac12b^*G^\dagger b. \tag{S11}
$$

Since \(\ker G=\ker D_Z\), \(b\perp\ker G\), so the minimum exists despite redundant coordinates. The gain applies to the specified full quadratic and paid space, without implying the cheapest construction or a nonlinear convergence rate.

## S5. Certificate, constraints, and model-error boundaries

For any local proposal \(s\), \(d=Hs+h=H(s-s_*)\). If \(q=d-Hz\), expand \(q^*H^{-1}q\) to obtain

$$
E_0^2=d^*H^{-1}d=2d^*z-z^*Hz+q^*H^{-1}q. \tag{S12}
$$

The original-proposal lower bound is \(\max(0,2d^*z-z^*Hz)\). The original upper adds \(\|q\|^2/\mu\), for verified \(H\succeq\mu I\). The repaired step \(s-z\) has only the tail error \(q^*H^{-1}q\). Adding the already removed energy to its acceptance test would evaluate the wrong step.

If \(M_0^{1/2}HM_0^{1/2}\succeq aI\), then \(q^*H^{-1}q\le q^*M_0q/a\). In the experiments, \(\Lambda=\lambda I\) gives \(\mu=\lambda\). The stored thin reduced factor \(Z_0\) satisfies \(H_0=\mathbb R_{\rm material}(Z_0^*Z_0)+\lambda I\); its Frobenius norm gives the conservative endpoint \(a=\lambda/(\lambda+\|Z_0\|_F^2)\). A Ritz value alone is not an enclosing endpoint. The unscaled indicator \(d^*M_0d\) remains an estimate.

For \(p=-\Phi(s-z)>0\), the optimal reduction is \(p_*=p+E_{\rm rem}^2/2\). Hence \(U_{\rm rem}\le2\tau p/(1-\tau)\) implies \(p/p_*\ge1-\tau\). The comparison uses \(\tau=0.01\). Floating-point certificates are conditional on the declared operators, endpoints, and action accuracy, not interval-arithmetic enclosures.

Constraints add a term absent from the unconstrained closure. If exact constrained minimizers satisfy \(H_is_i+h_i+n_i=0\), with \(n_i\) in the normal cone to the common feasible set, then

$$
H_0(s_C-s_0)=-\Delta Hs_C-D^*r-(n_C-n_0). \tag{S13}
$$

The normal-cone difference need not lie in the paired support. If both minimizers are known to lie on a common affine active face \(s=a+Fu\), one can apply the unconstrained algebra to the reduced coordinates with Jacobians \(J_iF\), residuals \(r+J_ia\), and the restricted regularizer, retaining any resulting residual-change terms. Active-face changes do not inherit that statement. No such constrained closure is claimed for the deployed projection routine.

The reused routine imposes \(\operatorname{Re}\chi\ge-0.5\) and \(\operatorname{Im}\chi\ge0\). Voxel coordinates are clipped. Gaussian coordinates use a Euclidean nearest-step QP in the physical tangent, solved by SLSQP with \(\texttt{ftol}=10^{-12}\) and 200 iterations. A violation above \(10^{-7}\) raises an error; the later line search has a \(10^{-8}\) feasibility allowance. This is not an \(H\)-metric constrained GN solve or a certified KKT solution. Projected proposals undergo full nonlinear Armijo acceptance; non-descent projections can trigger explicit projected-gradient fallback.

Likewise, replacing the numerical Maxwell operator or data changes \(H,h\) and the target minimizer. A certificate for one frozen discrete model says nothing by itself about that additional error. The finer-DDA check in S8 illustrates the distinction without providing a continuum error bound or an independent solver validation.

## S6. Existing campaign and numerical comparison contract

The experiment implementation was frozen before the held-out runs. It forms Jacobian-vector and adjoint actions rather than a dense full material Jacobian in the large cases. The six preflight identities were checked on a \(4^3\) DDA model with 192 complex current coordinates, six illuminations, 54 real material coordinates, base rank 16, and a two-complex-column bank. Recorded relative errors are \(3.881358\times10^{-15}\), \(8.329538\times10^{-14}\), \(6.540705\times10^{-14}\), \(4.270382\times10^{-16}\), \(1.260824\times10^{-16}\), and \(1.118307\times10^{-15}\) for T1–T6. They check implementation identities, not physical performance.

The main scene configurations are:

| Object | Grid | Tangent | Family | \(k\) | Max Re \(\chi\) | Max Im \(\chi\) | Illumination |
|---|---:|---|---|---:|---:|---:|---|
| 11 | 12 | Gaussian | safe | 1 | 0.490 | 0.063 | balanced |
| 12 | 12 | voxel | smooth | 2 | 1.826 | 0.109 | balanced |
| 13 | 12 | Gaussian | contact | 2.5 | 3.349 | 0.202 | balanced |
| 14 | 12 | Gaussian | partial | 2.2 | 1.793 | 0.142 | reduced diversity |
| 15 | 12 | Gaussian | curved | 2 | 1.211 | 0.084 | balanced |
| 16 | 14 | voxel | smooth | 2.1 | 1.233 | 0.097 | balanced |
| 17 | 14 | Gaussian | contact | 2.5 | 2.909 | 0.209 | balanced |
| 18 | 14 | voxel | partial | 2.3 | 1.257 | 0.086 | reduced diversity |

Contrast is \(\epsilon_r-1\). “Contact” denotes overlapping smooth inclusions, not a sharp internal dielectric junction. Balanced acquisition has three directions with two polarizations each; reduced diversity has six directions with one polarization each, not a limited aperture. All \(12^3\) inversions use \(14^3\) generated data; the three \(14^3\) stress objects use same-grid data. These are same-family discretization comparisons.

E2 uses unanchored Gs24 on the normalized base proposal and three fixed random directions. E3–E5 append full primal currents, with actual base rank 30. Residual candidates come from the dominant current correction directions of the base material proposal; Gs candidates use the next two discarded measurement modes; the TSOM-inspired bank uses two retained-orthogonalized directions from a bounded approximate domain bank. It is not the native TSOM optimizer. Generic material banks are Gaussian-random orthonormal columns. The E3 Krylov control starts from \(-M_0d_0\), applies \(M_0H\), and reorthogonalizes to the requested material rank; all 72 controls attain rank 48.

For the finite reference \(s_{\rm ref}\), let its remaining normal residual be \(d_{\rm ref}\) and \(\rho=\|d_{\rm ref}\|/\sqrt\lambda\). If a candidate's measured distance is \(e=\|s-s_{\rm ref}\|_H\), its true finite-model error is conditionally enclosed by

$$
[\max(0,e-\rho),e+\rho]. \tag{S14}
$$

Nonoverlapping intervals establish strict order. Relative ratios require both errors to be reference resolved; the near-floor Krylov errors therefore give bounded absolute orderings, not enormous measured gain factors. Figure 3 displays norm ratios; the statistical material-repair tables use squared-norm ratios. The eight objects, rather than the 24 anchors or 72 bank rows, are the bootstrap unit. Existing intervals use 10,000 object resamples and are descriptive for this small configuration set.

The registered numerical closure check divides the residual \(H\)-distance to the material span by the \(H\)-norm of \(s_C-s_0\). Its tolerance is \(\max(10^{-10},10\times\text{full-state solve residual})\). Five records fail, all with the TSOM-inspired bank at anchor 0:

| Object | Normalized distance | Absolute \(H\)-distance | Step-change \(H\)-norm | Threshold |
|---|---:|---:|---:|---:|
| 11 | 1.0939666e-8 | 2.2031698e-15 | 2.0139279e-7 | 1e-10 |
| 12 | 1.1782990e-9 | 2.7878889e-14 | 2.3660285e-5 | 1e-10 |
| 15 | 6.2132273e-10 | 1.0607294e-15 | 1.7072117e-6 | 1e-10 |
| 16 | 1.2359778e-9 | 3.4281067e-14 | 2.7735989e-5 | 1e-10 |
| 18 | 6.2408257e-10 | 3.1206828e-14 | 5.0004325e-5 | 1e-10 |

The other 67 pass. Small absolute distances do not waive the normalized failures. All paired-versus-direct optimization orderings can still hold, because they test a different quantity. The measured maximum Schur condition number, 1.061, is not evidence of a continuum distance from resonance.

## S7. Certificate calibration and complete-cost detail

The 96 original proposals consist of the base and three paired-bank proposals at each anchor. The saved calibration separates original-proposal safe recall from acceptance after executing the correction.

| Estimator | Safe recall | False safe | Coverage | Median \(H\)-norm effectivity | p90 | Maximum |
|---|---:|---:|---:|---:|---:|---:|
| Crude positive-curvature bound | 72.3% | 0 | 100.0% | 3.526 | 7.583 | 25.39 |
| Preconditioned indicator | 100.0% | 0 | 91.7% | 1.027 | 1.089 | 1.204 |
| Conservative spectral bound | 46.8% | 0 | 100.0% | 25.58 | 32.34 | 33.09 |
| Tail-corrected bound | 93.6% | 0 | 100.0% | 1.105 | 5.806 | 25.39 |

Effectivity is an upper error-norm estimate divided by the finite-reference proxy \(H\)-error norm; reference-induced effectivity intervals remain in the raw data. The tail method accepts 44 of 47 safe original proposals, rejects three safe proposals, and falsely accepts none of 49 unsafe proposals. All 96 proposals become acceptable after repair, a different event. Its mean correction count is 1.365 and mean recorded calibration time 0.5498 s. The unscaled preconditioned indicator has coverage only 88/96; zero observed false-safe decisions do not turn it into a bound.

![Existing certificate calibration](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A12/A12_FIGURES/Figure4_certificate.png)

*Fig. S1. Existing A12 Figure 4. Left: original squared errors and bounds. Right: upper \(H\)-norm effectivity and paid wall time. Zero lower bounds remain in the data but cannot appear on a logarithmic axis.*

The separate unpreconditioned eight-action Gauss/Radau diagnostic has 192 action rows. Its conditional bounds are supported by the eligible finite-reference intervals. Final \(H\)-norm effectivity is about 1.005, but the median paired time ratio to accepted B0 is 1.514. It is separately paid offline work, not an extra free consequence of the deployed PCG trajectory.

All E5 routes use the same \(M_0\), full actions, initial proposal, and 99% reduction target. B1 is a bounded material-history/Krylov adaptation; it compresses accumulated paid step/gradient history by SVD to eight directions and is not a reproduction of a named published solver. B2 is a selective sensitivity/adjoint adaptation, likewise not a full reproduction of adaptive-RB literature. B3 uses paired warm repair; B4 and B5 use paired and generic augmentation respectively. B4/B5 comparisons isolate bank source under the same route; B3/B4 also change execution route.

| Route | Method-stage ratio [95% CI] | Anchor-inclusive ratio [95% CI] | Cold-equivalent ratio [95% CI] | Median PCG corrections | Median full RHS |
|---|---:|---:|---:|---:|---:|
| B0 | 1 | 1 | 1 | 2 | 66 |
| B1 | 3.489 [3.457, 3.541] | 1.386 [1.248, 1.392] | 1.316 [1.187, 1.320] | 0 | 234 |
| B2 | 2.308 [2.264, 2.321] | 1.199 [1.123, 1.205] | 1.163 [1.093, 1.168] | 2 | 150 |
| B3 | 3.757 [3.528, 4.147] | 1.416 [1.330, 1.445] | 1.341 [1.235, 1.368] | 1 | 342 |
| B4 | 9.786 [9.430, 10.031] | 2.299 [1.908, 2.365] | 2.066 [1.678, 2.121] | 1 | 648 |
| B5 | 9.505 [9.287, 9.753] | 2.270 [1.880, 2.361] | 2.039 [1.657, 2.111] | 1 | 642 |

The B3 anchor-inclusive ratio is precisely 1.4162928067, with object-bootstrap interval [1.3295889298, 1.4454515137]. Each ratio first pairs the same object, anchor, and repeat with B0; six ratios per object are collapsed by a median, followed by an eight-object median and bootstrap. “Anchor-inclusive” is method wall time plus shared ROM setup and full-state anchor time. The cold equivalent additionally allocates common geometry setup across three anchors and must not be summed over methods as actual campaign elapsed time. Offline reference and truth generation are separate.

JVP/VJP call counts alone do not count batched columns. Full RHS counts each full-solve right-hand-side column. Candidate generation, adjoint pullbacks, full \(JZ/HZ\), QR, small solves, and refreshes remain charged. E3 timings omit candidate generation for direct steps and random generation/QR for generic banks, while single-sided tests are charged the shared paired build. Those unequal subtotals cannot support end-to-end claims.

Development measured a 1.9249 s median bank-plus-warm phase and an estimated PCG saving of −0.0762 s after removing that expense. The frozen controller disabled unpaid new banks; held-out timings were not used to learn a new rule. \(M_0\) is rebuilt at each anchor. No verified cross-anchor spectral transfer, dynamic rolling cost adaptation, or N-Port rank-saving result is part of the implemented evidence.

## S8. Nonlinear, noise, and finer-discretization evidence

Five preselected objects compare B0, B1, static ROM, and MURR. The other three provide B0 trajectories and anchors. All methods share \(\chi_0=0.1+0.04i\), prior weight \(10^{-5}\), initial LM damping 0.01 multiplied by 0.3 every three iterations, an 18-iteration cap, physical projection, and full nonlinear Armijo acceptance. “Static” denotes a fixed Gs24 mode-selection rule plus the current primal anchor, with operators refreshed at each state, not a single immutable operator.

The registered endpoint criteria are \(e_{\chi,m}\le1.02e_{\chi,B0}+0.002\) and \(e_{y,m}\le1.02e_{y,B0}+0.0002\). Receiver holdout uses a rotated receiver array at the same frequency and incident fields. This is distinct from object-level holdout from development. A B0 material error above 25% triggers an insufficient-recovery screen; passing it alone does not prove accurate or generally robust reconstruction.

| Object | B0 material error | Static material error | B0 seconds | Static seconds | B0/static state trials | Static endpoint noninferiority |
|---|---:|---:|---:|---:|---:|---|
| 11 | 2.810% | 3.906% | 94.7 | 73.4 | 18/18 | fail |
| 13 | 0.973% | 1.135% | 92.4 | 206.8 | 16/57 | pass |
| 15 | 0.447% | 0.540% | 101.8 | 140.9 | 18/38 | pass |
| 16 | 9.539% | 16.320% | 553.5 | 684.8 | 42/64 | fail |
| 17 | 0.067% | 0.128% | 233.8 | 371.2 | 16/34 | pass |

MURR follows B0 and has bitwise-identical saved materials in all five objects. Its nonlinear time ratio is 1.000 [0.994, 1.019]. This is a skip-path result, not paired enrichment evidence. Object 16's combined 9.54% B0/MURR material error contains real-part error 8.38% and imaginary-part error 36.82%; the static imaginary error is 50.33%. Visual loss and boundary errors are retained, not hidden behind the combined norm. A voxel tangent of dimension 5488 with 1536 real measurements has unregularized local Jacobian nullity at least 3952; that dimensional fact is not a global nonlinear nonuniqueness theorem.

![Existing reconstruction slices](/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A12/A12_FIGURES/Figure6_material_slices.png)

*Fig. S2. Existing A12 Figure 6. Truth, B0, static, and MURR on shared row scales. Object 13 uses finer-grid data; object 16 is same-grid stress. Slices are the positive central cell planes, not exactly \(z=0\). B0/MURR identity shows the skip route.*

The full nonlinear B0 trajectories reach as many as 37 PCG corrections as LM damping changes. The fixed-curvature E5 median cannot represent the entire nonlinear solve. Explicit gradient fallbacks, projection distances, and all accepted/rejected full-state trials remain paid.

Two Gaussian objects use one fixed complex noise realization at each nonzero level, shared across B0/MURR. The clean synthetic norm remains the normalization; these are controlled checks, not noise-blind deployment. The shared endpoints are:

| Object | Noise | Material error | Clean receiver-holdout relative error | B0 seconds | Linear certified proposals | Full-state trials | Gradient fallbacks |
|---|---:|---:|---:|---:|---:|---:|---:|
| 11 | 0% | 2.810% | 0.00005953 | 94.7 | 18/18 | 18 | 0 |
| 11 | 1% | 12.116% | 0.001587 | 596.5 | 18/18 | 164 | 11 |
| 11 | 3% | 18.879% | 0.004846 | 631.2 | 18/18 | 175 | 13 |
| 13 | 0% | 0.973% | 0.0007953 | 92.4 | 16/16 | 16 | 0 |
| 13 | 1% | 1.527% | 0.001801 | 94.5 | 16/16 | 16 | 0 |
| 13 | 3% | 2.362% | 0.004622 | 96.5 | 16/16 | 16 | 0 |

For object 11 at 1% noise, all 18 preprojection proposals were accepted by the local certificate, but only seven retain that provenance after 11 gradient replacements. At 3%, five retain it after 13 replacements. Neither count certifies constrained optimality. None of the four nonzero-noise endpoints crosses 25% material error, but skipped banks establish no paired-bank noise advantage.

The finer-grid endpoint check uses object 13 only. Saved Gaussian coefficients are lifted to \(14^3\) using the same continuous dictionary, with maximum fitting residual \(2.147\times10^{-15}\), followed by forward evaluation only. Exact B0/MURR material identity allows one shared fine solve.

| Method | Material error at \(12^3\) | Material error at \(14^3\) | Clean holdout at \(12^3\) | Clean holdout at \(14^3\) | Training field at \(14^3\) |
|---|---:|---:|---:|---:|---:|
| B0 / MURR | 0.009734 | 0.009733 | 0.0007953 | 0.007911 | 0.007910 |
| Static | 0.01135 | 0.01134 | 0.0008284 | 0.007947 | 0.007946 |

All entries are relative norms, not percentages. The 0.7809% coarse/fine truth-field discrepancy uses the training acquisition and fine truth-field norm. The approximately 0.00361 percentage-point static–B0 fine holdout difference uses the rotated receiver operator and its own truth norm. Their operators and denominators differ; the former is not an error bar or a resolution threshold for the latter. The evaluation reveals discretization dependence without providing a physical-error enclosure for method differences or a rigorous continuum bound from two grids. It was registered after initial held-out inspection, does not alter trajectories or thresholds, and remains within the same DDA family.

## S9. Evidence provenance and final scope

All large computations reported here had already completed on an RTX 4060 Laptop GPU, 8 GB, with an i9-14900HX host. They are historical campaign facts, not computation performed for this writing. The 20 main jobs and one fine-grid endpoint evaluation completed. Peak GPU allocation/reservation was 5.064/6.074 GiB. A separate development native-process failure with exit code 3221225477 and its successful same-configuration retry remain distinct records; the cause is not inferred from the exit code.

The immutable A12 experiment ledger has 1698 data rows, including 31 process records, 24 anchors, 35 nonlinear records, 576 iteration rows, 648 E3 method rows, 96 calibrations, and 288 E5 executions. It includes development and rejected/failed records; it is not a count of independent held-out experiments. The final writing read all rows and retained the original files.

Three analysis corrections remain essential. First, the legacy method-level constant \(\texttt{full_fallbacks}=0\) is not a count of projected-gradient fallback; per-iteration history is authoritative. Second, history objectives are pre-step while material errors are post-step; final objective is reconstructed from the last accepted trial or unchanged state rather than a rejected last attempt. Third, original linear certification and retained certificate provenance after projection/fallback are distinct events. These corrections affect reporting, not frozen algorithms or measurements.

The manuscript reuses original figure assets. Its Fig. 1/2/3/4 map to A12 Figure 1/2/3/5; supplementary Fig. S1/S2 map to A12 Figure 4/6. Corresponding CSV sidecars and per-anchor JSON remain authoritative. Raw figure ratios, empirical table values, finite-reference error intervals, complete histories, and source provenance are linked in the final claim ledger. Earlier conditional cross-anchor, TSOM-safe, and N-Port extensions remain in the historical supplement but are not added as contributions to this mechanism paper.

The manuscript and this supplement state and prove the theoretical results directly; numerical support is exclusively the existing campaign. The original S-DBIM/FFT-TSOM full-text gaps remain in the final reference audit. Inaccessible text is not evidence that an equivalent result does not exist.
