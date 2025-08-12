# Python Scripts for Class Vertical

## Equilibrium:

- *Equilibrium.py*:  Plots safety-factor, q, magnetic shear, s, higher order shear, s_2, pressure gradient p',
                     and derivate of pressure gradient, p'', versus radius, r
- *Profile.py*:	     Plots shaping functions P1, P1, Sigma, and S5 versus r				 

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
- *Hmat.py*:         Visualizes no-wall vacuum response matrix, H_mm'
- *iHmat.py*:        Visualizes inverse no-wall vacuum response matrix, iH_mm'

## Wall matrices:

- *Wall.py*:	     Visualizes wall matrices R_m^m and S_m^m
- *Imat.py*:	     Visualizes inverse of wall matrix I_m^m (should be Hermitian)
- *Gmat.py*:         Visualizes perfect-wall response matrix, G_mm'
- *iGmat.py*:        Visualizes inverse perfect-wall response matrix, iG_mm'
- *Bmat.py*:         Visualizes wall response matrix, B_mm'
- *Cmat.py*:         Visualizes wall response matrix, C_mm' (should be Hermtian)

## ODE Solution:

- *Matrix.py*:      Plots coupling matrices versus r
- *Solutions.py*:   Plots components of m-dominant solution vector versus r
- *Solution.py*:    Plots mth harmonic of components of m-dominant solution vector versus r
- *h.py*:           Plots adaptive ode integration data

## Ideal stability:

- *Umat.py*:		Visualizes total ideal energy matrix
- *Ideale.py*:      Plots poloidal harmonics of ideal eigenfunctions versus r
- *Ideale1.py*:     Plots kth poloidal harmonic of ideal eigenfunctions versus r
- *Evals.py*:	 	Plots eigenvalues of W, V, U matrices 
- *deltaW.py*:	 	Plots delta W values versus eigenfunction number
- *deltaW1.py*:	 	Plots select number of delta W values versus eigenfunction number
- *deltaW2.py*:	 	Plots delta W values versus eigenfunction number in range of y
- *yZSurface.py*:   Plots y and Z on plasma boundary associated with ideal eigenfunctions

## Visualization:

- *y.py*:           Plots y components of ideal eigenfunction in R, Z plane 
- *yZ.py*:          Plots y and Z components of ideal eignenfunction in R, Z plane
