# version 0.1.1 - 29/12/2021 this program plots the percentage of the XZ plane satisfying a limit on X
# it works with analize_data 0.8
# the files are hardcoded

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

datadir="pickled_data_final/"
infiles=(datadir+"results_Au_Elab_80_b6.pickle.gz",datadir+"results_Au_Elab_80_b6_10K.pickle.gz")
labels=("Au+Au, E$_{lab}$=80 AGeV, b=6fm, 1.08K ev.","Au+Au, E$_{lab}$=80 AGeV, b=6fm, 10.8K ev.")

colors=("blue","red","orange")
lines=("solid","dashed","dotted")

# parameters for the plots
plt.rcParams.update({'font.size': 16})
fig_size = plt.rcParams["figure.figsize"]
plt.rc('legend',fontsize=12)
fig_size[0] = 11
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size

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
    print('Syntax: ./make_plots_Au_Elab80_b6_comparison.py <Xlim> <output directory>')
    sys.exit(1)

Xlim=float(sys.argv[1])
od = sys.argv[2]
if (not os.path.exists(od)):
    os.mkdir(od)


area_percent=[]

for nif, ifile in enumerate(infiles):
    if (ifile[-3:] == ".gz"):
        print("Opening gzipped file " + ifile)
        fx = gzip.open(ifile, "rb")
    else:
        print("Opening file " + ifile)
        fx = open(ifile, "rb")

    data = pickle.load(fx)
    fx.close()

    lattice, tt, elim, Xlim_ebe, Ylim_ebe, XT_lim, RP_lim, N_events, en_per_part, Tmunu, X, Y, XT, RP, \
    N_cells_avg, N_cells_sd, Non_empty_cells_avg, Non_empty_cells_sd, N_cells_avg_fraction, \
    N_cells_sd_fraction,N_clumps_avg, N_clumps_percent_avg, N_clumps_single_cell_percent_tot, \
    N_clumps_single_cells_tot,N_clumps_empty_tot = data[:]
    
    nt=len(tt)

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_arr=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>pmin)

    if (nif==0):
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

    area_percent.append(np.zeros(nt,dtype=np.float64))

    for h in range(nt):
        area_good=0.
        area_system=0.
        for i in range(nx):
            for k in range(nz):
                if ((Tmunu[h,iT00,i,j0,k]>0) and (den[h,i,j0,k]>pmin)):
                    area_system=area_system+1
                    if (X_arr[h,i,j0,k]<Xlim):
                        area_good=area_good+1
        if (area_system>0):
            area_percent[nif][h]=area_good/area_system*100

sp="       "
tp='{:5.3f}'
vp='{:9.6f}'
xstring='{:3.2f}'.format(Xlim)

# percentage of the system lying on the reaction plane with X<0.3 
fout=open(od + "/X_lt_"+xstring+"_XZ_plane_percentage_Au_Elab80_b6.dat","w")
fout.write('# Au+Au at Elab = 80 AGeV, b = 6 fm y=z=0\n')
fout.write("# time [fm]      % area on y=0 plane (X < "+xstring+")\n")
plt.xlabel('t [fm]')
plt.ylabel("% area on y=0 plane (X < "+xstring+")")
plt.xlim(0,15)
plt.ylim(0,100)
for q in range(len(infiles)):
    plt.plot(tt, area_percent[q], color=colors[q], linewidth=2, linestyle=lines[q],label=labels[q])
    fout.write("\n\n# "+labels[q]+"\n\n")
    for i in range(len(tt)):
        fout.write(tp.format(tt[i])+sp+vp.format(area_percent[q][i])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(od + "/X_lt_"+xstring+"_XZ_plane_percentage_Au_Elab80_b6.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_lt_"+xstring+"_XZ_plane_percentage_Au_Elab80_b6.pdf", pad_inches=0.)
plt.close('all')
