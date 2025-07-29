import scipy.constants as sc

# Equilibrium
R_0 = 6.2
B_0 = 5.3
mu_0 = sc.mu_0

# 3, 1 surface
Ls    = 1.756
bh    = 0.0540
ab    = 2.588
ac    = 0.5816
Gboot = 6.381
Geccd = 1.433
Jmax  = bh * (ab - ac) * Gboot /Geccd
jmax  = (B_0 /mu_0 /R_0 /Ls) * Jmax

print ("3, 1 surface: Jmax = %11.4e jmax = %11.4e" % (Jmax, jmax))

# 2,1 surface
Ls    = 1.209
bh    = 0.0293
ab    = 3.061
ac    = 1.000
Gboot = 6.381
Geccd = 1.433
Jmax  = bh * (ab - ac) * Gboot /Geccd
jmax  = (B_0 /mu_0 /R_0 /Ls) * Jmax

print ("2, 1 surface: Jmax = %11.4e jmax = %11.4e" % (Jmax, jmax))
