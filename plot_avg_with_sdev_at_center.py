# production version 0.1.1 - 17/12/2021 this program plots the constraints vs time in the center of the grid and saves the plotted data in ascii files

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import matplotlib
import gzip
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm

N_input_files = len(sys.argv) - 1
if N_input_files != 1:
    print('Syntax: ./plot_avg_at_center_w_stdev.py <output directory>')
    sys.exit(1)

# parameters for the plots
plt.rcParams.update({'font.size': 16})
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 11
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
plt.rc('legend',fontsize=12)
colors = ("forestgreen","blue","red")
#linestyles = ('solid','dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1)))
linestyles = ('solid','dashed','dashdot')
colors_band = ("greenyellow","blueviolet","darkorange")

labels=["Ag+Ag, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.58AGeV, b=0-2.44fm ","Au+Au, "+r"$\mathrm{E}_{\mathrm{lab}}$"+"=1.23AGeV, b=0-3.3fm ","Au+Au, "+r"$\sqrt{s_{\mathrm{NN}}}$"+"=7.7GeV, b=0-3.3fm "]

# minimum pressure for considering a certain cell part of the system (in GeV/fm^3)
pmin = 1.e-4

#correspondences between the indexes of the 1D Tmunu array and the energy momentum rank 2 tensor
iT00=0
iT01=1
iT02=2
iT03=3
iT10=1
iT11=4
iT12=5
iT13=6
iT20=2
iT21=5
iT22=7
iT23=8
iT30=3
iT31=6
iT32=8
iT33=9


od = sys.argv[1]
if (not os.path.exists(od)):
    os.mkdir(od)

fx = open("time_evolution_comparison_w_stdev_at_center.pickle","rb")
data = pickle.load(fx)
fx.close()

times_arr,events_arr,Tmunu,Xebe,Yebe,Xavg,Yavg = data[:]

# plotted time range
time_start=0
time_end=25

# formats for the text data of the plots
sp="       "
tf='{:5.3f}'
pf='{:11.9e}'


# now we make the plots
# energy density
plt.xlabel('t [fm]')
plt.ylabel(r'$\varepsilon$ [GeV/fm$^3$]')
plt.xlim(time_start, time_end)
fout=open(od + "/energy_density_evolution_at_grid_center_w_errorbands.dat","w")
fout.write("# time [f]       energy density [GeV/fm^3]     standard deviation of energy density [GeV/fm^3]\n")
for i, tt in enumerate(times_arr):    
    tend=len(tt) # time dimension of Tmunu can be larger than len(tt)
    plt.plot(tt,Tmunu[i,0,:tend,iT00], color=colors[i], linewidth=2, linestyle=linestyles[i], label=labels[i])
    plt.fill_between(tt, Tmunu[i,0,:tend,iT00]-Tmunu[i,1,:tend,iT00], Tmunu[i,0,:tend,iT00]+Tmunu[i,1,:tend,iT00], alpha=0.2, 
             edgecolor=colors_band[i], facecolor=colors[i], linewidth=1, linestyle='dotted', antialiased=True)
    for mm in range(tend):
        if ((tt[mm]>=time_start) and (tt[mm]<=time_end)):
            val=Tmunu[i,0,mm,iT00]
            stdev=Tmunu[i,1,mm,iT00]
            fout.write(tf.format(tt[mm])+sp+pf.format(val)+sp+pf.format(stdev)+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/energy_density_evolution_at_grid_center_w_errorbands.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/energy_density_evolution_at_grid_center_w_errorbands.pdf", pad_inches=0.)
plt.close('all')



plt.xlabel('t [fm]')
plt.ylabel("Pressure anisotropy $X_{ebe}$")
plt.xlim(time_start, time_end)
plt.ylim(0,2)
fout=open(od + "/X_ebe_evolution_at_grid_center_w_errorbands.dat","w")
fout.write("# time [f]       X       standard deviation of X\n")
for i, tt in enumerate(times_arr):    
    tend=len(tt) # time dimension of Tmunu can be larger than len(tt)
    plt.plot(tt,Xebe[i,0,:tend], color=colors[i], linewidth=2, linestyle=linestyles[i], label=labels[i])
    plt.fill_between(tt, Xebe[i,0,:tend]-Xebe[i,1,:tend], Xebe[i,0,:tend]+Xebe[i,1,:tend], alpha=0.2, 
             edgecolor=colors_band[i], facecolor=colors[i], linewidth=1, linestyle='dotted', antialiased=True)
    for mm in range(tend):
        if ((tt[mm]>=time_start) and (tt[mm]<=time_end)):
            val=Xebe[i,0,mm]
            stdev=Xebe[i,1,mm]
            fout.write(tf.format(tt[mm])+sp+pf.format(val)+sp+pf.format(stdev)+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/X_ebe_evolution_at_grid_center_w_errorbands.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_ebe_evolution_at_grid_center_w_errorbands.pdf", pad_inches=0.)
plt.close('all')



plt.xlabel('t [fm]')
plt.ylabel("Off-diagonality $Y_{ebe}$")
plt.xlim(time_start, time_end)
plt.ylim(0,2)
fout=open(od + "/Y_ebe_evolution_at_grid_center_w_errorbands.dat","w")
fout.write("# time [f]       Y       standard deviation of Y\n")
for i, tt in enumerate(times_arr):    
    tend=len(tt) # time dimension of Tmunu can be larger than len(tt)
    plt.plot(tt,Yebe[i,0,:tend], color=colors[i], linewidth=2, linestyle=linestyles[i], label=labels[i])
    plt.fill_between(tt, Yebe[i,0,:tend]-Yebe[i,1,:tend], Yebe[i,0,:tend]+Yebe[i,1,:tend], alpha=0.2, 
             edgecolor=colors_band[i], facecolor=colors[i], linewidth=1, linestyle='dotted', antialiased=True)
    for mm in range(tend):
        if ((tt[mm]>=time_start) and (tt[mm]<=time_end)):
            val=Yebe[i,0,mm]
            stdev=Yebe[i,1,mm]
            fout.write(tf.format(tt[mm])+sp+pf.format(val)+sp+pf.format(stdev)+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/Y_ebe_evolution_at_grid_center_w_errorbands.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/Y_ebe_evolution_at_grid_center_w_errorbands.pdf", pad_inches=0.)
plt.close('all')

