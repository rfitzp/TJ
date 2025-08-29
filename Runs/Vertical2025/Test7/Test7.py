import matplotlib.pyplot as plt
import pandas as pd

df  = pd.read_csv("Results.out",  skiprows=0, delim_whitespace=True)
df1 = pd.read_csv("Results1.out", skiprows=0, delim_whitespace=True)
df2 = pd.read_csv("Results2.out", skiprows=0, delim_whitespace=True)
#df3 = pd.read_csv("Results3.out", skiprows=0, delim_whitespace=True)

fig = plt.figure (figsize = (12.0, 8.0))
plt.rc ('xtick', labelsize = 15) 
plt.rc ('ytick', labelsize = 15) 

plt.subplot (2, 2, 1)

plt.xlim (-1.05, 1.05)

plt.plot (df2.iloc[:,0], df2.iloc[:,2], color = 'red',    linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)=-0.5$')
plt.plot (df1.iloc[:,0], df1.iloc[:,2], color = 'green',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.0$')
plt.plot (df.iloc [:,0],  df.iloc[:,2], color = 'blue',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.5$')
plt.axvline (0.,                        color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$V_2(a)$',        fontsize = "15")
plt.ylabel (r'$\delta W_{nw}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 2)

plt.xlim (-1.05, 1.05)

plt.plot (df2.iloc[:,0], df2.iloc[:,3], color = 'red',    linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)=-0.5$')
plt.plot (df1.iloc[:,0], df1.iloc[:,3], color = 'green',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.0$')
plt.plot ( df.iloc[:,0],  df.iloc[:,3], color = 'blue',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.5$')
plt.axvline (0.,                        color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$V_2(a)$',        fontsize = "15")
plt.ylabel (r'$\delta W_{pw}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 3)

plt.xlim (-1.05, 1.05)

plt.plot (df2.iloc[:,0], df2.iloc[:,4], color = 'red',    linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)=-0.5$')
plt.plot (df1.iloc[:,0], df1.iloc[:,4], color = 'green',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.0$')
plt.plot ( df.iloc[:,0],  df.iloc[:,4], color = 'blue',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.5$')
plt.axvline (0.,                        color = 'black',  linewidth = 1.5, linestyle = 'dotted')

plt.xlabel (r'$V_2(a)$',   fontsize = "15")
plt.ylabel (r'$\alpha_w$', fontsize = "15")
plt.legend (fontsize = "12")

plt.subplot (2, 2, 4)

plt.xlim (-1.05, 1.05)

plt.plot (df2.iloc[:,0], df2.iloc[:,5], color = 'red',    linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)=-0.5$')
plt.plot (df1.iloc[:,0], df1.iloc[:,5], color = 'green',  linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.0$')
plt.plot ( df.iloc[:,0],  df.iloc[:,5], color = 'blue',   linewidth = 1,   linestyle = 'dotted', marker = 's', fillstyle = 'none', markersize = 5, label = r'$H_3(a)= 0.5$')
plt.axvline (0.,                        color = 'black',  linewidth = 1.5, linestyle = 'dotted')
 
plt.xlabel (r'$V_2(a)$',       fontsize = "15")
plt.ylabel (r'$\hat{\gamma}$', fontsize = "15")
plt.legend (fontsize = "12")

plt.tight_layout ()

#plt.show ()    
plt.savefig ("Fig8.pdf")
