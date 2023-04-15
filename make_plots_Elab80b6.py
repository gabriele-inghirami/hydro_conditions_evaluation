# version 0.1.1 - 18/12/2021 this program plots the constraints vs time in the center of the grid
# it works with analize_data 0.8

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

# minimum value of the sum of the diagonal term to consider a point part of the system [GeV/fm^3]
pmin=0.0001

# parameters for the plots
plt.rcParams.update({'font.size': 16})
fig_size = plt.rcParams["figure.figsize"]
plt.rc('legend',fontsize=12)
fig_size[0] = 11
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
colors = ("crimson", "darkorange", "forestgreen")
linestyles=("solid","dashed","dotted")

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

# we get the name output directory
N_input_files = len(sys.argv) - 1

if N_input_files != 2:
    print('Syntax: ./make_plots_for_article.py <input pickled file> <output directory>')
    sys.exit(1)

ifile = sys.argv[1]
od = sys.argv[2]
if (not os.path.exists(od)):
    os.mkdir(od)

if (ifile[-3:] == ".gz"):
    print("Opening gzipped file " + ifile)
    fx = gzip.open(ifile, "rb")
else:
    print("Opening file " + ifile)
    fx = open(ifile, "rb")

data = pickle.load(fx)
fx.close()

lattice, tt, elim, Xlim, Ylim, XT_lim, RP_lim, N_events, en_per_part, Tmunu, X, Y, XT, RP, \
N_cells_avg, N_cells_sd, Non_empty_cells_avg, Non_empty_cells_sd, N_cells_avg_fraction, \
N_cells_sd_fraction,N_clumps_avg, N_clumps_percent_avg, N_clumps_single_cell_percent_tot, \
N_clumps_single_cells_tot,N_clumps_empty_tot = data[:]
    
nt=len(tt)

den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
X_arr=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>pmin)

num=3*abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:])
Y_arr=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>pmin)

coord=[]
for i in range(0,3):
    coord.append(np.linspace(lattice["origin"][i]+lattice["spacing"][i]/2., lattice["origin"][i]+(lattice["dimensions"][i]-1)*lattice["spacing"][i]+lattice["spacing"][i]/2.,num=lattice["dimensions"][i], endpoint=True))
xx,yy,zz=np.array(coord[:])


nx=len(xx)
ny=len(yy)
nz=len(zz)

i0=np.argmin(abs(xx-0))
j0=np.argmin(abs(yy-0))
k0=np.argmin(abs(zz-0))

i2=np.argmin(abs(xx-2))
i6=np.argmin(abs(xx-6))

points=[i0,i2,i6]
labels=["x=y=z=0","x=2,y=z=0","x=6,y=z=0"]

area_percent=np.zeros(nt,dtype=np.float64)

for h in range(nt):
    area_good=0.
    area_system=0.
    for i in range(nx):
        for k in range(nz):
            if ((Tmunu[h,iT00,i,j0,k]>0) and (den[h,i,j0,k]>pmin)):
                area_system=area_system+1
                if (X_arr[h,i,j0,k]<0.3):
                    area_good=area_good+1
    if (area_system>0):
        area_percent[h]=area_good/area_system*100

sp="       "
tp='{:5.3f}'
vp='{:9.6f}'
# X (avg Tmunu)
fout=open(od + "/X_evolution_Au_Elab80_b6.dat","w")
fout.write('# Au+Au at Elab = 80 AGeV, b = 6 fm, 1080 events, y=z=0\n')
fout.write("# time [fm]        X (x=0)       X (x=2)      X (x=6)\n\n")
plt.xlabel('t [fm]')
plt.ylabel("Pressure anisotropy $X$")
plt.ylim(0,2)
plt.xlim(0,15)
for i, p in enumerate(points): 
    plt.plot(tt, X_arr[:,p,j0,k0], color=colors[i], linewidth=2, linestyle=linestyles[i], label=labels[i])
for i in range(len(tt)):
    fout.write(tp.format(tt[i]))
    for p in range(len(points)):
        fout.write(sp+vp.format(X_arr[i,points[p],j0,k0]))
    fout.write("\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#999999', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/X_evolution_Au_Elab80_b6.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_evolution_Au_Elab80_b6.pdf", pad_inches=0.)
plt.close('all')

# Y (avg Tmunu)
fout=open(od + "/Y_evolution_Au_Elab80_b6.dat","w")
fout.write('# Au+Au at Elab = 80 AGeV, b = 6 fm, 1080 events, y=z=0\n')
fout.write("# time [fm]        Y (x=0)       Y (x=2)      Y (x=6)\n\n")
plt.xlabel('t [fm]')
plt.ylabel("Off-diagonality $Y$")
plt.xlim(0,15)
for i, p in enumerate(points): 
    plt.plot(tt, Y_arr[:,p,j0,k0], color=colors[i], linewidth=2, linestyle=linestyles[i], label=labels[i])
for i in range(len(tt)):
    fout.write(tp.format(tt[i]))
    for p in range(len(points)):
        fout.write(sp+vp.format(Y_arr[i,points[p],j0,k0]))
    fout.write("\n")
fout.close()
plt.ylim(0,0.3)
plt.minorticks_on()
plt.grid(b=True, color='#999999', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/Y_evolution_Au_Elab80_b6.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/Y_evolution_Au_Elab80_b6.pdf", pad_inches=0.)
plt.close('all')

# percentage of the system lying on the reaction plane with X<0.3 
fout=open(od + "/X_lt_03_XZ_plane_percentage_Au_Elab80_b6.dat","w")
fout.write('# Au+Au at Elab = 80 AGeV, b = 6 fm, 1080 events, y=z=0\n')
fout.write("# time [fm]      % area on y=0 plane (X < 0.3)\n\n")
plt.xlabel('t [fm]')
plt.ylabel("% area on y=0 plane (X < 0.3)")
plt.xlim(0,15)
plt.ylim(0,100)
plt.plot(tt, area_percent, color="blue", linewidth=2, linestyle="solid",label="Au+Au, E$_{lab}$=80AGeV, b=6fm")
for i in range(len(tt)):
    fout.write(tp.format(tt[i])+sp+vp.format(area_percent[i])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#999999', linestyle='-')
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(od + "/X_lt_03_XZ_plane_percentage_Au_Elab80_b6.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_lt_03_XZ_plane_percentage_Au_Elab80_b6.pdf", pad_inches=0.)
plt.close('all')
