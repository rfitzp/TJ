# Python Scripts for Class TJ

## Vacuum matrices:

- *Metric.py*:       Plots metric data at plasma boundary
- *Vacuum.py*:       Visualizes vacuum matrices P_m^m' and R_m^m
- *Vacuum1.py*:      Visualizes vacuum matrices Q_m^m' and S_m^m
- *PRmat.py*:        Visualizes vacuum matrix PR^mm' (should be Hermitian)
- *QSmat.py*:        Visualizes vacuum matrix QS^mm' (should be Hermitian)
- *PSmat.py*:        Visualizes vacuum matrix PS^mm' (should be unit matrix)
- *QPmat.py*:        Visualizes vacuum matrix QP^mm' (should be Hermitian)
- *RSmat.py*:        Visualizes vacuum matrix RS^mm' (should be Hermitian)
- *SPmat.py*:        Visualizes vacuum matrix SP^mm' (should be unit matrix)
- *Hmat.py*:         Visualizes vacuum response matrix, H_mm'
- *iHmat.py*:        Visualizes inverse vacuum response matrix, iH_mm'

## Wall matrices:

- *Wall.py*:	     Visualizes wall matrices R_m^m and S_m^m
- *Imat.py*:	     Visualizes inverse of wall matrix I_m^m (should be Hermitian)
- *Gmat.py*:         Visualizes wall response matrix, G_mm'
- *iGmat.py*:        Visualizes inverse wall response matrix, iG_mm'

## ODE Solution:

- *Matrix.py*:      Plots coupling matrices versus r
= *km.py*:          Plots km values versus r
- *Solutions.py*:   Plots components of m-dominant solution vector versus r
- *Solution.py*:    Plots mth harmonic of components of m-dominant solution vector versus r
- *h.py*:           Plots adaptive ode integration data

## Tearing eigenfunctions:

- *Full.py*:         Plots poloidal harmonics of fully reconnected solution vector versus r
- *Unrc.py*:         Plots poloidal harmonics of unreconnected solution vector versus r
- *Unrc1.py*:        Plots kth poloidal harmonic of unreconnected solution vector versus r
- *Psi.py*:          Plots psi components of unreconnected solution vector in R, Z plane  
- *PsiZ.py*:         Plots psi and Z components of unreconnected solution vector in R, Z plane

## Electromagnetic torques:

- *Tfull.py*:       Plots angular momentum flux associated with pair of fully reconnected solutions versus r
- *Tunrc.py*:       Plots angular momentum flux associated with pair of unreconnected solutions versus r

## Tearing dispersion relation:

- *Emat.py*:	 Visualizes tearing stability matrix

## RMP coils:

- *Psix.py*:           Plots poloidal harmonics of RMP at plasma boundary
- *Chi.py*:            Plots RMP drives at rational surfaces
- *Rmp.py*:            Plots ideal response of plasma to RMP
- *Rmp1.py*:           Plots kth poloidal harmonics of ideal response of plasm to RMP
- *Psir.py*:           Plots psi components of ideal response of plasma to RMP in R, Z plane
- *PsiZr.py*:          Plots psi and Z components of ideal response of plasma to RMP in R, Z plane
- *PsiRmpSurface.py*:  Plots RMP data on plasma boundary
- *Gamma.py*:          Plots expansions of Psi_x and Psi^rmp at plasma boundary in terms of ideal eigenfunctions

## Ideal stability:

- *Ideal.py*:		Plots poloidal harmonics of ideal solutions launched from magnetic axis versus r
- *Ideal1.py*:		Plots kth poloidal harmomic of ideal solutions launched from magnetic axis versus r
- *Internal.py*:        Plots poloidal harmonics of internal ideal solutions launched from magnetic axis versus r
- *Umat.py*:		Visualizes total ideal energy matrix
- *Ideale.py*:      	Plots polidal harmonics of ideal eigenfunctions versus r
- *Ideale1.py*:     	Plots kth poloidal harmonic of ideal eigenfunctions versus r
- *Evals.py*:	 	Plots eigenvalues of W, V, U matrices 
- *deltaW.py*:	 	Plots delta W values versus eigenfunction number
- *Xi.py*:              Plots Xi Fourier harmonics on plasma boundary
- *PsiJSurface.py*: 	Plots psi and J on plasma boundary associated with ideal eigenfunctions
- *lambda.py*:          Plots eigenvalues of plasma energy matrix versus r
- *crit.py*:		Plots smallest eigenvalue of plasma energy matrix versus r