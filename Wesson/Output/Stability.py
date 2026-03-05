import matplotlib.pyplot as plt
import math
import netCDF4 as nc

infile = open("Stability.out", "r") 

qa  = []
nu1 = []
nu2 = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qa .append(float(c1))
    nu1.append(float(c2))
    nu2.append(float(c3))
 
fig = plt.figure(figsize=(10.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)

plt.xlim(qa[0], qa[-1])
plt.ylim(1., max(nu1))
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qa, nu1, color='black',  linewidth=2, linestyle = 'dashed')
plt.plot(qa, nu2, color='black',  linewidth=2, linestyle = 'dashed')

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
#plt.legend(fontsize="20")

plt.show()
#plt.savefig("width.eps")
    
