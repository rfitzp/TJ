import matplotlib.pyplot as plt
import math

infile = open("m2n1r12.out", "r") 

q1 = []
r1 = []
d1 = []
w1 = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    c4 = numbers[3]
    q1.append(float(c1))
    r1.append(float(c2))
    d1.append(float(c3))
    w1.append(float(c4))


r1.insert(0, 1.)
w1.insert(0, 0.)
d1.insert(0, 0.)
    
infile = open("m2n1r11.out", "r") 

q2 = []
r2 = []
d2 = []
w2 = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    c4 = numbers[3]
    q2.append(float(c1))
    r2.append(float(c2))
    d2.append(float(c3))
    w2.append(float(c4))

infile = open("m2n1r10.out", "r") 

q3 = []
r3 = []
d3 = []
w3 = []

for line in infile: 

    numbers = line.split() 
    c1 = numbers[0] 
    c2 = numbers[1] 
    c3 = numbers[2] 
    c4 = numbers[3]
    q3.append(float(c1))
    r3.append(float(c2))
    d3.append(float(c3))
    w3.append(float(c4))      

fig = plt.figure(figsize=(10.0, 8.0))
plt.rc('xtick', labelsize=17) 
plt.rc('ytick', labelsize=17) 

plt.subplot(1, 1, 1)

plt.xlim(0.3, 1.01)
plt.ylim(0., 0.5)
#plt.axvline(1.0, color='black', linewidth=1., linestyle='dotted')

plt.plot(r3, w3, color='black',  linewidth=2, linestyle = 'solid',  label=r'$\bar{b}=1.0$')
plt.plot(r2, w2, color='black',  linewidth=2, linestyle = 'dashed', label=r'$\bar{b}=1.1$')
plt.plot(r1, w1, color='black',  linewidth=2, linestyle = 'dotted', label=r'$\bar{b}=1.2$')

plt.ylabel(r'$\overline{W}_s$', fontsize="20")
plt.xlabel(r'$\bar{r}_s$',      fontsize='20')
plt.legend(fontsize="20")

#plt.show()

plt.savefig("Figure11_4.pdf")
    
