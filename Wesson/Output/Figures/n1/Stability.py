import matplotlib.pyplot as plt
import math
import netCDF4 as nc

fna = 'm2rw100.nc'
dsa = nc.Dataset(fna)
qaa = dsa['q_xa']
nua = dsa['nu']

fnb = 'm2rw11.nc'
dsb = nc.Dataset(fnb)
qab = dsb['q_xa']
nub = dsb['nu']

fnc = 'm3rw100.nc'
dsc = nc.Dataset(fnc)
qac = dsc['q_xa']
nuc = dsc['nu']

fnd = 'm3rw11.nc'
dsd = nc.Dataset(fnd)
qad = dsd['q_xa']
nud = dsd['nu']

fne = 'm4rw100.nc'
dse = nc.Dataset(fne)
qae = dse['q_xa']
nue = dse['nu']

fnf = 'm4rw11.nc'
dsf = nc.Dataset(fnf)
qaf = dsf['q_xa']
nuf = dsf['nu']

infile = open("wm2rw100.out", "r") 

qax  = []
nu1x = []
nu2x = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qax .append(float(c1))
    nu1x.append(float(c2))
    nu2x.append(float(c3))
 
infile = open("wm2rw11.out", "r") 

qay  = []
nu1y = []
nu2y = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qay .append(float(c1))
    nu1y.append(float(c2))
    nu2y.append(float(c3))

infile = open("wm3rw100.out", "r") 

qaz  = []
nu1z = []
nu2z = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qaz .append(float(c1))
    nu1z.append(float(c2))
    nu2z.append(float(c3))
 
infile = open("wm3rw11.out", "r") 

qau  = []
nu1u = []
nu2u = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qau .append(float(c1))
    nu1u.append(float(c2))
    nu2u.append(float(c3))

infile = open("wm4rw100.out", "r")     

qav  = []
nu1v = []
nu2v = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qav .append(float(c1))
    nu1v.append(float(c2))
    nu2v.append(float(c3))   

infile = open("wm4rw11.out", "r")     

qaw  = []
nu1w = []
nu2w = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    qaw .append(float(c1))
    nu1w.append(float(c2))
    nu2w.append(float(c3))       
    
fig = plt.figure(figsize=(10.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(3, 2, 1)

plt.xlim(1., 8.)
plt.ylim(1., 6.)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qax,  nu1x, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qax,  nu2x, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qaa,  nua,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=2, n=1~~~\bar{b}=10$", fontsize = 15)
plt.yticks([1, 3, 5])
#plt.legend(fontsize="20")

plt.subplot(3, 2, 2)

plt.xlim(1., 8.)
plt.ylim(1., 6.)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qay,  nu1y, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qay,  nu2y, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qab,  nub,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=2, n=1~~~\bar{b}=1.1$", fontsize = 15)
plt.yticks([1, 3, 5])
#plt.legend(fontsize="20")

plt.subplot(3, 2, 3)

plt.xlim(2., 12.)
plt.ylim(1., 4.)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qaz,  nu1z, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qaz,  nu2z, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qac,  nuc,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xticks([3, 6, 9, 12])

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=3, n=1~~~\bar{b}=10$", fontsize = 15)
#plt.legend(fontsize="20")

plt.subplot(3, 2, 4)

plt.xlim(2., 12.)
plt.ylim(1., 4.)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qau,  nu1u, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qau,  nu2u, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qad,  nud,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xticks([3, 6, 9, 12])

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=3, n=1~~~\bar{b}=1.1$", fontsize = 15)
#plt.legend(fontsize="20")

plt.subplot(3, 2, 5)

plt.xlim(3.5, 4.5)
plt.ylim(1., 1.5)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qav,  nu1v, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qav,  nu2v, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qae,  nue,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xticks([3.5, 4., 4.5])

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=4, n=1~~~\bar{b}=10$", fontsize = 15)
#plt.legend(fontsize="20")

plt.subplot(3, 2, 6)

plt.xlim(3.5, 4.5)
plt.ylim(1., 1.5)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(qaw,  nu1w, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qaw,  nu2w, color = 'black', linewidth = 2, linestyle = 'dashed')
plt.plot(qaf,  nuf,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xticks([3.5, 4., 4.5])

plt.ylabel(r'$\nu$',   fontsize="20")
plt.xlabel(r'$q_a$',   fontsize='20')
plt.title(r"$m=4, n=1~~~\bar{b}=1.1$", fontsize = 15)
#plt.legend(fontsize="20")

plt.tight_layout ()

#plt.show()

plt.savefig("Figure11_2.pdf")
    
