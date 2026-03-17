import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

df  = pd.read_csv("Results.out",  skiprows=0, delim_whitespace=True)
df1 = pd.read_csv("Results1.out", skiprows=0, delim_whitespace=True)
df2 = pd.read_csv("Results2.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 1, 1)

plt.xlim (1.0, 1.5)

plt.plot (df.iloc [:,0],  df.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'solid',  label = r'$p_0=0.05$')
plt.plot (df1.iloc[:,0], df1.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'dashed', label = r'$p_0=0.10$')
plt.plot (df2.iloc[:,0], df2.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'dotted', label = r'$p_0=0.20$')

plt.xlabel (r'$\bar{b}$',      fontsize = "20")
plt.ylabel (r'$\hat{\gamma}$', fontsize = "20")
plt.legend (fontsize = "15")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure12_7.pdf")
