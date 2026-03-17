import matplotlib.pyplot as plt
import pandas as pd

df  = pd.read_csv("Results.out",  skiprows=0, delim_whitespace=True)
df1 = pd.read_csv("Results1.out", skiprows=0, delim_whitespace=True)
df2 = pd.read_csv("Results2.out", skiprows=0, delim_whitespace=True)
df3 = pd.read_csv("Results3.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (12.0, 6.0))
plt.rc ('xtick', labelsize = 20) 
plt.rc ('ytick', labelsize = 20) 

plt.subplot (1, 2, 1)

plt.xlim (0.0, 0.21)

plt.plot (df1.iloc[:,0], df1.iloc[:,2], color = 'black',  linewidth = 2,   linestyle = 'solid',   label = r'$E_a=0.5$')
plt.plot (df.iloc [:,0],  df.iloc[:,2], color = 'black',  linewidth = 2,   linestyle = 'dashed',  label = r'$E_a=1.0$')
plt.plot (df2.iloc[:,0], df2.iloc[:,2], color = 'black',  linewidth = 2,   linestyle = 'dotted',  label = r'$E_a=1.5$')
plt.plot (df3.iloc[:,0], df3.iloc[:,2], color = 'black',  linewidth = 2,   linestyle = 'dashdot', label = r'$E_a=2.0$')

plt.xlabel (r'$\epsilon_a$',          fontsize = "20")
plt.ylabel (r'$\delta \hat{W}_{nw}$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (1, 2, 2)

plt.xlim (0.0, 0.21)

plt.plot (df1.iloc[:,0], df1.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'solid',   label = r'$E_a=0.5$')
plt.plot (df.iloc [:,0],  df.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'dashed',  label = r'$E_a=1.0$')
plt.plot (df2.iloc[:,0], df2.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'dotted',  label = r'$E_a=1.5$')
plt.plot (df3.iloc[:,0], df3.iloc[:,5], color = 'black',  linewidth = 2,   linestyle = 'dashdot', label = r'$E_a=2.0$')

plt.xlabel (r'$\epsilon_a$',   fontsize = "20")
plt.ylabel (r'$\hat{\gamma}$', fontsize = "20")
plt.legend (fontsize = "15")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Figure12_5.pdf")
