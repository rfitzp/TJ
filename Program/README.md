# Program TJ

Program to calculate tearing stability matrix and tearing eigenfunctions for an aspect-ratio
expanded toroidal plasma equilibrium. The program deals with a perturbed magnetic field consisting
of a single toroidal harmonic and a range of different coupled poloidal harmonics. The program
solves the equations of marginally-stable ideal magnetohydrodynamics thoughout the plasma. These equations
become singular at various rational magnetic flux-surfaces within the plasma. The solution is
integrated to just before each rational surface and then jumped across the surface using an analytic
solution that is valid in the immediate vicinity of the surface. Either free or fixed boundary conditions
are imposed at the plasma boundary. The free boundary conditions are generated via an expansion in toroidal
functions in the vacuum region surrounding the plasma. 