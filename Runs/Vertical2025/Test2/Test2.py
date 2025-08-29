import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("Results.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (-2.05, 2.05)

plt.plot (df.iloc[:,0], df.iloc[:,2], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$m_{max}=5$')
plt.plot (df.iloc[:,0], df.iloc[:,7], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = '$m_{max}=10$')
plt.axhline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2(a)$',        fontsize = "15")
plt.ylabel (r'$\delta W_{nw}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 2)

plt.xlim (-2.05, 2.05)

plt.plot (df.iloc[:,0], df.iloc[:,3], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = '$m_{max}=5$')
plt.plot (df.iloc[:,0], df.iloc[:,8], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = '$m_{max}=10$')
plt.axvline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2(a)$',        fontsize = "15")
plt.ylabel (r'$\delta W_{pw}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 3)

plt.xlim (-2.05, 2.05)

plt.plot (df.iloc[:,0], df.iloc[:,4], color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = r'$m_{max}=5$')
plt.plot (df.iloc[:,0], df.iloc[:,9], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = r'$m_{max}=10$')
plt.axvline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2(a)$',   fontsize = "15")
plt.ylabel (r'$\alpha_w$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 4)

plt.xlim (-2.05, 2.05)

plt.plot (df.iloc[:,0], df.iloc[:,5],  color = 'blue',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 10, label = r'$m_{max}=5$')
plt.plot (df.iloc[:,0], df.iloc[:,10], color = 'red',   linewidth = 1,   linestyle = 'dotted', marker = '^', fillstyle = 'none', markersize = 10, label = r'$m_{max}=10$')
plt.axhline (0.,                       color = 'black', linewidth = 1.5, linestyle = 'dotted')
plt.axvline (0.,                      color = 'black', linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$H_2(a)$',       fontsize = "15")
plt.ylabel (r'$\hat{\gamma}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.tight_layout ()

#plt.show ()
plt.savefig ("Fig1.pdf")

