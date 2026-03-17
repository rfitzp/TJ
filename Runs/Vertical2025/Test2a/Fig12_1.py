import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("Results.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 2, 1)

plt.xlim (0.0, 2.0)

plt.plot (df.iloc[:,0], df.iloc[:,7], color = 'black',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                      color = 'black',   linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$E_a$',                 fontsize = "20")
plt.ylabel (r'$\delta \hat{W}_{nw}$', fontsize = "20")

plt.subplot (1, 2, 2)

plt.xlim (0.0, 2.0)

plt.plot (df.iloc[:,0], df.iloc[:,10], color = 'black',   linewidth = 2,   linestyle = 'solid')
plt.axhline (0.,                       color = 'black',   linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$E_a$',          fontsize = "20")
plt.ylabel (r'$\hat{\gamma}$', fontsize = "20")

plt.tight_layout ()

#plt.show ()
plt.savefig ("Figure12_3.pdf")

