# Python Scripts for Program TJ


## Plasma equilibrium:

	- Flux.py        - Plots r, theta coordinate system in R, Z plane
	- Fluxw.py       - Plots r, omega coordinate system in R, Z plane
	- Equilibrium.py - Plots components of aspect-ratio expanded equilibrium versus r
	- Safety.py      - Plots safety-factor, q, magnetic shear, s, higher order shear, s_2, pressure gradient and
	  		   derivate of pressure gradient, p'', versus radius, r
	- Shear.py       - Plots magnetic shear, amd related functions, versus r
	- Shape.py       - Plots values of shaping functions at plasma boundary
	- Shaping.py     - Plots shaping functions versus r
	- Profile.py     - Plots shaping function, S_1, and profiles function, P_1, P_2, P_3, versus radius, r
	- Boundary.py    - Plots data relating to plasma boundary

## Vacuum matrices:

   	- Vacuum.py      - Visualizes vacuum matrices P_m^m' and R_m^m
	- Amat.py        - Visualizes vacuum matrix A_m^m'
	- Hmat.py        - Visualizes vacuum response matrix, H_m^m'
	- Metric.py      - Plots metric data at plasma boundary

## ODE Solution:

	- Matrix.py      - Plots L_m^m', M_m^m', N_m^m', and P_m^m' coupling matrices versus r
	- Solutions.py   - Plots psi and Z components of m-dominant solution vector versus r
	- Solution.py    - Plots mth harmonic of psi and Z components of m-dominant solution vector versus r
	- h.py           - Plots adaptive ode integration step-length, h, and truncation error, err, versus radius, r

## Tearing eigenfunctions:

	- Full.py        - Plots poloidal harmonics of psi and Z Fourier components of fully reconnected solution vector associated with
	  		   given rational surface versus r
	- Unrc.py        - Plots poloidal harmonics of psi and Z Fourier components of unreconnected solution vector associated with
	  		   given rational surface versus r
	- Unrc1.py       - Plots kth poloidal harmonic of psi and Z components of unreconnected solution vector associated with given rational surface versus r.
	- Psi.py         - Plots psi components of unreconnected solution vector associated with given rational surface in R, Z plane
	- PsiZ.py        - Plots psi and Z components of unreconnected solution vector associated with given rational surface in R, Z plane

## Electromagnetic torques:

	- Tfull.py       - Plots angular momentum flux associated with pair of fully reconnected solutions versus r
	- Tunrc.py       - Plots angular momentum flux associated with pair of unreconnected solutions versus r

