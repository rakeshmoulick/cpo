import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
from celluloid import Camera
indx = 45000
fig = plt.figure(figsize=(7,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

file_name1 = 'ec%d.txt'%(indx)
file_name2 = 'eh%d.txt'%(indx)
#file_name3 = 'i%d.txt'%(indx)
path = '../data/'

# Load cold electron phase space
if os.path.exists(pjoin(path,file_name1)):
    xc, yc = np.loadtxt(pjoin(path,file_name1),unpack=True)
else:
    print('No Data')
    exit()

# Load hot electrons phase space
if os.path.exists(pjoin(path,file_name2)):
    xh, yh = np.loadtxt(pjoin(path,file_name2),unpack=True)
else:
    print('No Data')
    exit()

		
plt.scatter(xh, yh, s=0.3, alpha=0.5, color='r')
plt.scatter(xc, yc, s=0.3, alpha=0.5, color='b')	

#plt.xlim([-0.0005, 0.0025])
#plt.ylim([-3E6, 3E6])
plt.xlabel('x')
plt.ylabel('v')
#plt.title(indx)
#plt.grid(True)	
path_fig = './'
plt.savefig(pjoin(path_fig,'part.png'))
plt.show()