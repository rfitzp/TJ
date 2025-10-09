import numpy as np
import matplotlib.pyplot as plt

file_names='./slayer/HRii_Vz.out';
with open(file_names, 'r') as f:
    header=f.readline().strip().split();
    print(header)

    slayer=[];
    for line in f:
        row = list(map(float, line.strip().split()));
        slayer.append(row);

file_named='./data/HRii_data.txt';
with open(file_named, 'r') as f:
    header=f.readline().strip().split();
    print(header)

    data=[];
    for line in f:
        row = list(map(float, line.strip().split()));
        data.append(row);

slayer_array = np.array(slayer)
nslayer = len(slayer_array[:,0])

Qs = np.zeros(nslayer, dtype='float')
ReDs = np.zeros(nslayer, dtype='float')
ImDs = np.zeros(nslayer, dtype='float')

for i in range(nslayer):
    Qs[i] = slayer_array[i, 0];
    ReDs[i] = slayer_array[i, 2];
    ImDs[i] = slayer_array[i, 3];

data_array = np.array(data)
ndata = len(data_array[:,0])
print(np.shape(data_array))

Qd = np.zeros(ndata, dtype='float')
ReDd = np.zeros(ndata, dtype='float')
ImDd = np.zeros(ndata, dtype='float')

for i in range(ndata):
    Qd[i] = data_array[i, 0];
    ReDd[i] = data_array[i, 3];
    ImDd[i] = data_array[i, 4];

fig, ax = plt.subplots(1,2);

ax[0].plot(Qs, ReDs, label='slayer')
ax[0].plot(Qd, ReDd, '--', label='ps solver')
ax[0].set_xlabel('Q')
ax[0].set_ylabel(r'$\Re(\hat{\Delta})$')
ax[0].grid()
ax[1].plot(Qs, ImDs, label='slayer')
ax[1].plot(Qd, ImDd, '--', label='ps solver')
ax[1].set_xlabel('Q')
ax[1].set_ylabel(r'$\Im(\hat{\Delta})$')
ax[1].legend()
ax[1].grid()
plt.suptitle('HRii')
plt.tight_layout()
plt.show()
