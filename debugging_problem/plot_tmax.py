import numpy as np
import matplotlib.pyplot as plt
import os
import gmsh
import csv
from mpl_toolkits.axes_grid1 import host_subplot


plt.rcParams['text.usetex'] = True


main_dir = os.path.dirname(__file__)
os.chdir(main_dir)

tmax_abs =  np.genfromtxt('error/tmax_abs.txt', delimiter=' ', ndmin=2)
time = tmax_abs[:, 1]
tmax_abs = tmax_abs[:, 7]
tmax_rel =  np.genfromtxt('error/tmax_rel.txt', delimiter=' ', ndmin=2)[:,7]

fig, ax = plt.subplots()

ax.plot(time, tmax_abs, color='#440154', label='Absolute error')
ax.plot(time, tmax_rel, color='#21918c', linestyle='--',label='Relative error')

plt.xlim(0, 2)
# plt.ylim(1e-13, 1)
#ax_b.yaxis.get_major_locator().set_params(numticks=999)
#ax_b.yaxis.get_minor_locator().set_params(numticks=999, subs=[2, 4, 6, 8])

ax.set_xlabel("Time [s]")
ax.set_ylabel('Error in max($T$)[K]')
ax.legend(loc='best')
ax.grid(True, which='both', lw=0.5)
#plt.savefig("h-phi-b-err.pdf")
plt.show()