import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("Results.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (0.0, 0.21)

plt.plot (df.iloc[:,0], df.iloc[:,2], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$\delta W_{nw}(5)$')
plt.plot (df.iloc[:,0], df.iloc[:,7], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = '$\delta W_{nw}(10)$')
plt.axhline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0.0, 0.21)

plt.plot (df.iloc[:,0], df.iloc[:,3], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$\delta W_{pw}(5)$')
plt.plot (df.iloc[:,0], df.iloc[:,8], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = '$\delta W_{pw}(10)$')
plt.axhline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0.0, 0.21)

plt.plot (df.iloc[:,0], df.iloc[:,4], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = r'$\alpha_{w\,nw}(5)$')
plt.plot (df.iloc[:,0], df.iloc[:,9], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = r'$\alpha_{w\,pw}(10)$')
plt.axhline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon$', fontsize = "15")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0.0, 0.21)

plt.plot (df.iloc[:,0], df.iloc[:,5],  color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = r'$\hat{\gamma}_{w\,nw}(5)$')
plt.plot (df.iloc[:,0], df.iloc[:,10], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = r'$\hat{\gamma}_{w\,pw}(10)$')
plt.axhline (0.,                       color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$\epsilon$', fontsize = "15")
plt.legend (fontsize = "15")

plt.tight_layout ()

plt.show ()    

