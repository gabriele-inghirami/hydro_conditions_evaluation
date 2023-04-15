# version 1.5.1 - 27/12/2021 this program plots the percentage of cells of the fireball with respect to time, under different max values of
# pressure anisotropy and T^munu off diagonality

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
from matplotlib.pyplot import cm
#from itertools import cycle

matplotlib.use('Agg')

# minimum pressure for considering a certain point part of the system (in GeV/fm^3)
pmin = 1.e-4

# we get the name output directory
N_input_files = len(sys.argv) - 1

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_1D_plots_percentage_XY_vs_time.py <pickled input file list> <output directory>')
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
plt.rcParams.update({'legend.fontsize': 12})
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
colors_line = ("forestgreen","blue","red","orange","darkviolet")
linestyles = ('solid','dashed','dashdot','dotted',(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1)))
# we need color bands only for plots with three curves
colors_band = ("greenyellow","blueviolet","darkorange")

od = sys.argv[2]
if (not os.path.exists(od)):
    os.mkdir(od)

# correspondences between the indexes of the 1D Tmunu array and the energy momentum rank 2 tensor
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

# dictionary with the results in the case of average T^munu only
res_avgT={}
time_dict={}
max_nt=0

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

    n_eps = len(elim)
    n_X = len(Xlim)
    n_Y = len(Ylim)
    n_XT = len(XT_lim)
    n_RP = len(RP_lim)
    nt = len(tt)
    nx, ny, nz = lattice["dimensions"][:]

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    num=3.*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)

    taken_cells=np.zeros(nt,dtype=np.float64)
    syst_cells=np.zeros(nt,dtype=np.float64)
    
    plot_data=np.zeros((nt,n_X,n_Y,n_eps), dtype=np.float64)

    # now we make the plots
    # event by event case
    for x in range(n_X):
        Xval = '{:3.2f}'.format(Xlim[x])
        for y in range(n_Y):
                    Yval = '{:3.2f}'.format(Ylim[y])
                    plt.xlabel('t [fm]')
                    plt.ylabel('Fraction of the fireball')
                    plt.xlim(-3, tt[-1])
                    
                    fp=open(od + "/percentage_of_fireball_above_min_energy_density_vs_time_ebe_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                    fp.write("# X_ebe < "+Xval+", Y_ebe < "+Yval+"\n\n")

                    for i in range(n_eps):
                       epsval = '{:03d}'.format(int(elim[i] * 1000))
                       fp.write("# energy density > "+epsval+" MeV/fm^3\n")
                       fp.write("# time       volume fraction      std_dev volume fraction\n")
                       plt.plot(tt, N_cells_avg_fraction[:, i, x, y, 0, 0], color=colors_line[i], linewidth=2, linestyle=linestyles[i], label=r"$\varepsilon$ > " + epsval + "$MeV/fm^3$")
                       plt.fill_between(tt, N_cells_avg_fraction[:, i, x, y, 0, 0]-N_cells_sd_fraction[:, i, x, y, 0, 0], N_cells_avg_fraction[:, i, x, y, 0, 0]+N_cells_sd_fraction[:, i, x, y, 0, 0],
                         alpha=0.2, edgecolor=colors_band[i], facecolor=colors_line[i], linewidth=1, linestyle='dashdot', antialiased=True)
                       for it in range(nt):
                           fp.write('{:5.2f}'.format(tt[it])+"      "+'{:9.6e}'.format(N_cells_avg_fraction[it, i, x, y, 0, 0])+"      "+'{:9.6e}'.format(N_cells_sd_fraction[it, i, x, y, 0, 0])+"\n")
                       fp.write("\n\n")

                    plt.minorticks_on()
                    plt.grid(b=True, color='#BBBBBB', linestyle='-')
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_ebe_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".png", dpi=300, pad_inches=0.)
                    plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_ebe_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".pdf", pad_inches=0.)
                    plt.close('all')
                    fp.close()

    # average Tmunu case
    for x in range(n_X):
        Xval = '{:3.2f}'.format(Xlim[x])
        for y in range(n_Y):
            Yval = '{:3.2f}'.format(Ylim[y])
            plt.xlabel('t [fm]')
            plt.ylabel('Fraction of the fireball')
            plt.xlim(-3, tt[-1])
            fp=open(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
            fp.write("# X < "+Xval+", Y < "+Yval+" (from average T^munu)\n\n")
            for i in range(n_eps):
                taken_cells[:]=0.
                syst_cells[:]=0.
                epsval = '{:03d}'.format(int(elim[i] * 1000))

                fp.write("# energy density > "+epsval+" MeV/fm^3\n")
                fp.write("# time       volume fraction\n")

                for it in range(nt):
                    for ix in range(nx):
                        for iy in range(ny):
                            for iz in range(nz):
                                if (Tmunu[it,iT11,ix,iy,iz] + Tmunu[it,iT22,ix,iy,iz] + Tmunu[it,iT33,ix,iy,iz] < pmin):
                                    continue
                                syst_cells[it] = syst_cells[it] + 1
                                if ((X_avg[it,ix,iy,iz]<Xlim[x]) and (Y_avg[it,ix,iy,iz]<Ylim[y]) and (Tmunu[it,iT00,ix,iy,iz] >= elim[i])):
                                    taken_cells[it] = taken_cells[it] + 1

                perc=np.divide(taken_cells,syst_cells,out=np.zeros(taken_cells.shape, dtype=np.float64), where=syst_cells!=0)
                epsval = '{:03d}'.format(int(elim[i] * 1000))
                plt.plot(tt, perc, color=colors_line[i], linewidth=2, linestyle=linestyles[i], label=r"$\varepsilon$ > " + epsval + "$MeV/fm^3$")
                plot_data[:,x,y,i]=perc[:]
                for it in range(nt):
                    fp.write('{:5.2f}'.format(tt[it])+"      "+'{:9.6e}'.format(perc[it])+"\n")
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='upper right')
            plt.tight_layout()
            plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_" + key + "_Xlim_" + Xval + "_Ylim_" + Yval + ".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

    res_avgT[key]=plot_data
    time_dict[key]=tt
    if len(tt)>max_nt:
        max_nt=len(tt)
        tt_big=tt


# comparison plots between different systems
for x in range(n_X):
    Xval = '{:3.2f}'.format(Xlim[x])
    for y in range(n_Y):
        Yval = '{:3.2f}'.format(Ylim[y])
        for i in range(n_eps):
            epsval = '{:03d}'.format(int(elim[i] * 1000))
            fp=open(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.xlabel('t [fm]')
            plt.ylabel('Fraction of the fireball')
            plt.xlim(-3, tt[-1])
            fp.write("# X < "+Xval+", Y < "+Yval+", edens >= "+epsval+" MeV/fm^3 (from average T^munu)\n")
            fp.write("# column 1: time [fm]\n")
            line_number=2
            for key, value in infiles.items():
                fp.write("# column "+'{:2d}'.format(line_number)+": volume fraction of the system in "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            l=0
            for key, value in infiles.items():
                #plt.plot(time_dict[key], res_avgT[key][:,x,y,i], color=next(colors), linewidth=2, linestyle=next(linestyles), label=value[1])
                plt.plot(time_dict[key], res_avgT[key][:,x,y,i], color=colors_line[l], linewidth=2, linestyle=linestyles[l], label=value[1])
                l=l+1
            for it in range(max_nt):
                fp.write('{:5.2f}'.format(tt_big[it])+"      ")
                for key, value in infiles.items():
                    if it < len(time_dict[key]):
                        fp.write('{:9.6e}'.format(res_avgT[key][it,x,y,i])+"      ")
                    else:
                        fp.write('{:9.6e}'.format(0.)+"      ")
                fp.write("\n")
            plt.minorticks_on()
            #plt.legend(loc='lower right')
            plt.legend(loc='best')
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.tight_layout()
            plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/percentage_of_fireball_above_min_energy_density_vs_time_average_Tmunu_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()
