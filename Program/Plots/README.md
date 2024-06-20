# Python scripts for program TJ


## Equilibrium:

	- Flux.py        - Plots equilibrium magnetic flux-surfaces in R, Z plane
	- Equilibrium.py - Plots components of aspect-ratio expanded equilibrium versus r
	- Safety.py      - Plots safety-factor, q, magnetic shear, s, higher order shear, s_2, pressure gradient and derivate of pressure gradient, p'', versus radius, r
	- Shear.py       - Plots magnetic shear, amd related functions, versus r
	- Shape.py       - Plots values of shaping functions at plasma boundary
	- Shaping.py     - Plots shaping functions versus r
	- Profile.py     - Plots shaping function, S_1, and profiles function, P_1, P_2, P_3, versus radius, r

## Matrices:

	- Matrix.py      - Plots L_m^m', M_m^m', N_m^m', and P_m^m' coupling matrices versus r
	- VacRes.py      - Visualizes vacuum solution residual matrices A_m^m', B_m^m', and C_m^m'
	- Hmat.py        - Visualizes vacuum matching matrix, H_m^m'

## ODESolution:

	- h.py           - Plots adaptive ode integration step-length, h, and truncation error, err, versus radius, r
	- Solutions.py   - Plots psi and Z components of m-dominant solution vector versus r
	- Solution.py    - Plots mth harmonic of psi and Z components of m;-dominant solution vector versus r

## EigenFunctions:

	- Full.py        - Plots psi and Z components of fully reconnected solution vector associated with given rational surface versus r
        - Unrc.py        - Plots psi and Z components of unreconnected solution vector associated with given rational surface versus r
	- Unrc1.py       - Plots kth harmonic of psi and Z components of unreconnected solution vector associated with given rational surface versus r.

## Torques:

	- Tfull.py       - Plots angular momentum flux associated with pair of fully reconnected solutions versus r
        - Tunrc.py       - Plots angular momentum flux associated with pair of unreconnected solutions versus r