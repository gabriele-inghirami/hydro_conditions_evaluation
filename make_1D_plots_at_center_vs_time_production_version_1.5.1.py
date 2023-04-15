# production version 1.5.1 - 14/12/2021 this program plots the constraints vs time in the center of the grid and saves the plotted data in ascii files

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import matplotlib
import gzip
import gc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm

matplotlib.use('Agg')

# we get the name output directory
N_input_files = len(sys.argv) - 1

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_1D_plots_at_center_vs_time.py <pickled list of input files> <output directory>')
    sys.exit(1)

# path of the file containing the pickled list of input files
datafiles_list = sys.argv[1]

# directory of the file containing the pickled list of input files
datapath = os.path.split(datafiles_list)[0]+"/"

# we check if the directory with the data exists
if not os.path.exists(datafiles_list):
    print("Error, the file " + datafiles_list + " that should contain the list of input data files does not exist!")
    sys.exit(1)

# we use a dictionary for the files to be analized
# the key will enter in the output filenames
# the values are a tuple with: 0) the name of the files 1) a string with the common part of the plot titles
# this dictionary is save in a pickled datafile that it is used also by other python scripts
with open(datafiles_list,"rb") as inputfiles:
    infiles=pickle.load(inputfiles)

# parameters for the plots
plt.rcParams.update({'font.size': 16})
plt.rc('legend',fontsize=12)
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 11
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
colors = ("crimson", "darkorange", "olive", "forestgreen", "blue")
linestyles = ((0, ()),(0, (5, 5)),(0, (1, 1)),(0, (3, 1, 1, 1)),(0, (3, 1, 1, 1, 1, 1)))

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


od = sys.argv[2]
if (not os.path.exists(od)):
    os.mkdir(od)

eps_dict={}
en_per_part_dict={}
X_dict={}
Y_dict={}
X_ae_dict={}
Y_ae_dict={}
XT_dict={}
RP_dict={}
XT_ae_dict={}
RP_ae_dict={}
tt_dict={}

for key, value in infiles.items():
    if (value[0][-3:] == ".gz"):
        print("Opening gzipped file " + datapath + value[0])
        fx = gzip.open(datapath + value[0], "rb")
    else:
        print("Opening file " + datapath + value[0])
        fx = open(datapath + value[0], "rb")

    data = pickle.load(fx)
    fx.close()

    lattice, tt, elim, Xlim, Ylim, XT_lim, RP_lim, N_events, en_per_part, Tmunu, X, Y, XT, RP, \
    N_cells_avg, N_cells_sd, Non_empty_cells_avg, Non_empty_cells_sd, N_cells_avg_fraction, \
    N_cells_sd_fraction,N_clumps_avg, N_clumps_percent_avg, N_clumps_single_cell_percent_tot, \
    N_clumps_single_cells_tot,N_clumps_empty_tot = data[:]
    
    tt_dict[key]=np.copy(tt)

    i0=int(math.floor(lattice["dimensions"][0]/2))
    j0=int(math.floor(lattice["dimensions"][1]/2))
    k0=int(math.floor(lattice["dimensions"][2]/2))

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_ae_arr=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    X_ae_dict[key]=np.copy(X_ae_arr[:,i0,j0,k0])

    num=3*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_ae_arr=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    Y_ae_dict[key]=np.copy(Y_ae_arr[:,i0,j0,k0])

    eps_dict[key]=np.copy(Tmunu[:,iT00,i0,j0,k0])
    en_per_part_dict[key]=np.copy(en_per_part[:,i0,j0,k0])
    X_dict[key]=np.copy(X[:,i0,j0,k0])
    Y_dict[key]=np.copy(Y[:,i0,j0,k0])
    XT_dict[key]=np.copy(XT[:,i0,j0,k0])
    RP_dict[key]=np.copy(RP[:,i0,j0,k0])
   
    del lattice
    del tt
    del Xlim
    del Ylim
    del XT_lim
    del RP_lim
    del N_events
    del Tmunu
    del en_per_part
    del X
    del Y
    del XT
    del RP
    del N_cells_avg
    del N_cells_sd
    del Non_empty_cells_avg
    del Non_empty_cells_sd
    del N_cells_avg_fraction
    del N_cells_sd_fraction
    del N_clumps_avg
    del N_clumps_percent_avg
    del N_clumps_single_cell_percent_tot
    del N_clumps_single_cells_tot
    del N_clumps_empty_tot
    del data
    del num
    del Y_ae_arr
    del X_ae_arr
    del den
    gc.collect()

# plotted time range
time_start=0
time_end=25


# now we make the plots
# energy density

fout=open(od + "/energy_density_evolution_at_grid_center.dat","w")
sp="    "
#plt.title("Energy density in the center of the grid")
plt.xlabel('t [fm]')
plt.ylabel(r'$\varepsilon$ [GeV/fm$^3$]')
plt.xlim(time_start, time_end)
i=0
for key, value in infiles.items():    
    plt.plot(tt_dict[key], eps_dict[key], color=colors[i], linewidth=2, linestyle=linestyles[i], label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-energy density [GeV/fm^3]\n\n")
    for mm in range(len(tt_dict[key][:])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:11.9e}'.format(eps_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/energy_density_evolution_at_grid_center.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/energy_density_evolution_at_grid_center.pdf", pad_inches=0.)
plt.close('all')


# energy per baryon
#plt.title("Energy per net baryon in the center of the grid")
fout=open(od + "/energy_per_net_baryon_evolution_at_grid_center.dat","w")
plt.xlabel('t [fm]')
plt.ylabel(r'$\varepsilon$ [MeV]')
plt.xlim(time_start, time_end)
i=0
for key, value in infiles.items():    
    plt.plot(tt_dict[key], 1000*en_per_part_dict[key], color=colors[i], linewidth=2, linestyle=linestyles[i], label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-energy per net baryon [MeV]\n\n")
    for mm in range(len(tt_dict[key][:])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:10.6e}'.format(1000*en_per_part_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/energy_per_net_baryon_evolution_at_grid_center.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/energy_per_net_baryon_evolution_at_grid_center.pdf", pad_inches=0.)
plt.close('all')

# X (avg single events)
#plt.title("Pressure anisotropy (average of single events)")
fout=open(od + "/X_ebe_evolution_at_grid_center_avg_single_events.dat","w")
plt.xlabel('t [fm]')
plt.ylabel("Pressure anisotropy $X_{ebe}$")
plt.xlim(time_start, time_end)
plt.ylim(0,2)
i=0
for key, value in infiles.items():    
    if (key == "C_Elab_2_b_0_1_2"):
        time_stop = np.argmin(abs(np.array(tt_dict[key])-10))+1
    elif (key == "ArKCl_Elab_1_756_b_0_1_784"):
        time_stop = np.argmin(abs(np.array(tt_dict[key])-17))+1
    else:
        time_stop = np.argmin(abs(np.array(tt_dict[key])-time_end))+1

    plt.plot(tt_dict[key][:time_stop], X_dict[key][:time_stop], color=colors[i], linewidth=2, linestyle=linestyles[i], label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-X_ebe (average many single events)\n\n")
    for mm in range(len(tt_dict[key][:time_stop])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:10.6e}'.format(X_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/X_ebe_evolution_at_grid_center_avg_single_events.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_ebe_evolution_at_grid_center_avg_single_events.pdf", pad_inches=0.)
plt.close('all')

# Y (avg single events)
#plt.title("Off diagonality (average of single events)")
fout=open(od + "/Y_ebe_evolution_at_grid_center_avg_single_events.dat","w")
plt.xlabel('t [fm]')
plt.ylabel("Off-diagonality $Y_{ebe}$")
plt.xlim(time_start, time_end)
plt.ylim(0,1.75)
i=0
for key, value in infiles.items():    
    if (key == "C_Elab_2_b_0_1_2"):
        time_stop = np.argmin(abs(np.array(tt_dict[key])-10))+1
    elif (key == "ArKCl_Elab_1_756_b_0_1_784"):
        time_stop = np.argmin(abs(np.array(tt_dict[key])-17))+1
    else:
        time_stop = np.argmin(abs(np.array(tt_dict[key])-time_end))+1
    plt.plot(tt_dict[key][:time_stop], Y_dict[key][:time_stop], color=colors[i], linewidth=2, linestyle=linestyles[i], label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-Y_ebe (average many single events)\n\n")
    for mm in range(len(tt_dict[key][:time_stop])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:10.6e}'.format(Y_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/Y_ebe_evolution_at_grid_center_avg_single_events.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/Y_ebe_evolution_at_grid_center_avg_single_events.pdf", pad_inches=0.)
plt.close('all')

# X (avg Tmunu)
#plt.title("Pressure anisotropy (from average "+r"$T^{\mu\nu}$"+")")
fout=open(od + "/X_evolution_at_grid_center_avg_Tmunu.dat","w")
plt.xlabel('t [fm]')
plt.ylabel("Pressure anisotropy $X$")
plt.xlim(time_start, time_end)
plt.ylim(0,2)
i=0
for key, value in infiles.items():    
    plt.plot(tt_dict[key][:], X_ae_dict[key][:], color=colors[i], linewidth=2, linestyle=linestyles[i], label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-X (average T^{\mu\nu})\n\n")
    for mm in range(len(tt_dict[key][:])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:10.6e}'.format(X_ae_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/X_evolution_at_grid_center_avg_Tmunu.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/X_evolution_at_grid_center_avg_Tmunu.pdf", pad_inches=0.)
plt.close('all')

# Y (avg Tmunu)
#plt.title("Off diagonality (from average "+r"$T^{\mu\nu}$"+")")
fout=open(od + "/Y_evolution_at_grid_center_avg_Tmunu.dat","w")
plt.xlabel('t [fm]')
plt.ylabel("Off-diagonality $Y$")
plt.xlim(time_start, time_end)
plt.ylim(0,0.5)
i=0
for key, value in infiles.items():    
    plt.plot(tt_dict[key][:], Y_ae_dict[key][:], color=colors[i], linewidth=2, linestyle=linestyles[i],
                         label=value[1])
    i=i+1
    if (i>1):
         fout.write("\n\n\n")
    fout.write("# "+value[1]+"\n")
    fout.write("# 1-time [fm]      2-Y (average T^{\mu\nu})\n\n")
    for mm in range(len(tt_dict[key][:])):
        if ((tt_dict[key][mm]>=time_start) and (tt_dict[key][mm]<=time_end)):
            fout.write('{:5.3f}'.format(tt_dict[key][mm])+sp+'{:10.6e}'.format(Y_ae_dict[key][mm])+"\n")
fout.close()
plt.minorticks_on()
plt.grid(b=True, color='#BBBBBB', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig(od + "/Y_evolution_at_grid_center_avg_Tmunu.png", dpi=300, pad_inches=0.)
plt.savefig(od + "/Y_evolution_at_grid_center_avg_Tmunu.pdf", pad_inches=0.)
plt.close('all')
