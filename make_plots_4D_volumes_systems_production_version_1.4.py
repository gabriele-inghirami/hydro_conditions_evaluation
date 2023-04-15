# version 1.4 - 10/12/2021 this program plots the integrals of the results,
#               comparing all systems together for a fixed set of constraints

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import ScalarFormatter
from matplotlib.pyplot import cm
from itertools import cycle

# minimum pressure density in GeV/fm^3
pmin = 0.0001

# we get the name output directory
N_input_files = len(sys.argv) - 1

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_plots_4D_volumes.py <pickled list of input files> <output directory>')
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
fig_size[0] = 11
fig_size[1] = 6
plt.rcParams["figure.figsize"] = fig_size
#colors = ("crimson", "darkorange", "olive", "forestgreen", "blue", "deepskyblue")

def format_label(num, pos):
    decimal_digits=0
    if num != 0 :
        exponent = int(math.floor(math.log10(abs(num))))
        mant = round(num / float(10**exponent), decimal_digits)
        formatted_string = r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(mant, exponent, decimal_digits)
    else:
        formatted_string = "$0$"

    return formatted_string



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

od = sys.argv[2]
if (not os.path.exists(od)):
    os.mkdir(od)

# dictionary with the results 
res_ebe={}
res_avg={}
res_ebe_fraction_integral={}
res_avg_fraction_integral={}
res_avg_fraction_tot_integral={}
res_ebe_fraction_tot_integral={}
time_dict={}

first_time=True
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

    if first_time:
        nx,ny,nz = lattice["dimensions"][:]
        n_eps = len(elim)
        n_X = len(Xlim)
        n_Y = len(Ylim)

    plot_data_ebe=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)
    plot_data_ebe_percent=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)
    plot_data_ebe_fraction_integrals=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)
    plot_data_avg=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)
    plot_data_avg_percent=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)
    plot_data_fraction_integrals=np.zeros((n_X,n_Y,n_eps),dtype=np.float64)

    dt=tt[1]-tt[0]
    nt=len(tt)
    total_time=tt[-1]-tt[0]

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    num=3.*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
  
    for x in range(n_X):
    #for x in range(1,2):
      for y in range(n_Y):
      #for y in range(1,2):
        for i in range(n_eps):
        #for i in range(0,1):
          plot_data_ebe[x,y,i]=np.sum(N_cells_avg[:,i,x,y,0,0])*dt
          plot_data_ebe_percent[x,y,i]=np.sum(N_cells_avg_fraction[:,i,x,y,0,0])*dt
          denom=np.sum(Non_empty_cells_avg[:])
          if (denom>0):
              plot_data_ebe_fraction_integrals[x,y,i]=np.sum(N_cells_avg[:])/denom
          counts=0
          fraction_integral=0
          all_times_system_counts=0
          for it in range(nt):
               system_counts=0
               good_counts=0
               for ix in range(nx):
                   for iy in range(ny):
                       for iz in range(nz):
                           if (den[it,ix,iy,iz]>=pmin):
                               system_counts=system_counts+1
                               all_times_system_counts=all_times_system_counts+1
                               if((Tmunu[it,iT00,ix,iy,iz]>=elim[i]) and (X_avg[it,ix,iy,iz]<Xlim[x]) and (Y_avg[it,ix,iy,iz]<Ylim[y])):
                                   counts=counts+1
                                   good_counts=good_counts+1
               if (system_counts!=0):
                   fraction_integral=fraction_integral+good_counts/system_counts
          plot_data_avg[x,y,i]=counts*dt
          plot_data_avg_percent[x,y,i]=fraction_integral*dt
          if(all_times_system_counts!=0):
               plot_data_fraction_integrals[x,y,i]=counts/all_times_system_counts
     
    res_ebe[key]=plot_data_ebe
    res_avg[key]=plot_data_avg
    res_ebe_fraction_integral[key]=plot_data_ebe_percent/total_time
    res_avg_fraction_integral[key]=plot_data_avg_percent/total_time
    res_avg_fraction_tot_integral[key]=plot_data_fraction_integrals
    res_ebe_fraction_tot_integral[key]=plot_data_ebe_fraction_integrals

colors=cycle(cm.rainbow(np.linspace(0,1,len(res_ebe))))
mmm=cycle(['o','v','^','<','>','p','*',"P","X","D"])
marker_size=10
for x in range(n_X):
#for x in range(1,2):
    Xval = '{:3.2f}'.format(Xlim[x])
    for y in range(n_Y):
    #for y in range(1,2):
        Yval = '{:3.2f}'.format(Ylim[y])
        for i in range(n_eps):
        #for i in range(0,1):
            epsval = '{:03d}'.format(int(elim[i]*1000.))
            fp=open(od + "/4Dvolume_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('<4D volume> [fm$^4$]')
            fp.write("# X_ebe < "+Xval+", Y_ebe < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_ebe[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.xticks([])
            for key, value in infiles.items():
                plt.plot(line_number, res_ebe[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

            fp=open(od + "/4Dvolume_average_fraction_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('Average volume fraction [fm$^{-1}$]')
            fp.write("# Time averaged ratio between the volume satifying a set of constraints and the total volume of the system (i.e. with p>1e-4 GeV/fm^3)\n")
            fp.write("# X_ebe < "+Xval+", Y_ebe < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_ebe_fraction_integral[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.ylim(0,1)
            plt.gca().yaxis.set_major_formatter(ScalarFormatter())
            plt.xticks([])
            for key, value in infiles.items():
                plt.plot(line_number, res_ebe_fraction_integral[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_average_fraction_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_average_fraction_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

            fp=open(od + "/4Dvolume_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('<4D volume> [fm$^4$]')
            fp.write("# X < "+Xval+", Y < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_avg[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.gca().yaxis.set_major_formatter(FuncFormatter(format_label))
            plt.xticks([])
            for key, value in infiles.items():
                plt.plot(line_number, res_avg[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

            fp=open(od + "/4Dvolume_fraction_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('Average volume fraction [fm$^{-1}$]')
            fp.write("# Time averaged ratio between the volume satifying a set of constraints and the total volume of the system (i.e. with p>1e-4 GeV/fm^3)\n")
            fp.write("# X_ebe < "+Xval+", Y_ebe < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_avg_fraction_integral[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.ylim(0,1)
            plt.xticks([])
            plt.gca().yaxis.set_major_formatter(ScalarFormatter())
            for key, value in infiles.items():
                plt.plot(line_number, res_avg_fraction_integral[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='lower right')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_fraction_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_fraction_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

            fp=open(od + "/4Dvolume_ratio_of_integrals_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('( 4D volume constrained )/( 4D volume system )')
            fp.write("# Ratio between the time integrated volume satifying a set of constraints and the time integrated total volume of the system (i.e. with p>1e-4 GeV/fm^3)\n")
            fp.write("# X_ebe < "+Xval+", Y_ebe < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_ebe_fraction_tot_integral[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.ylim(0,1)
            plt.gca().yaxis.set_major_formatter(ScalarFormatter())
            plt.xticks([])
            for key, value in infiles.items():
                plt.plot(line_number, res_ebe_fraction_tot_integral[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='best')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_ratio_of_integrals_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_ratio_of_integrals_system_comparison_ebe_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()

            fp=open(od + "/4Dvolume_ratio_of_integrals_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".dat","w")
            plt.ylabel('( 4D volume constrained )/( 4D volume system )')
            fp.write("# Ratio between the time integrated volume satifying a set of constraints and the time integrated total volume of the system (i.e. with p>1e-4 GeV/fm^3)\n")
            fp.write("# X < "+Xval+", Y < "+Yval+", edens >= "+epsval+" MeV/fm^3\n")
            line_number=1
            for key, value in infiles.items():
                fp.write("# "+'{:2d}'.format(line_number)+"      "+value[1]+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            for key, value in infiles.items():
                fp.write('{:2d}'.format(line_number)+"  "+'{:12.2e}'.format(res_avg_fraction_tot_integral[key][x,y,i])+"\n")
                line_number=line_number+1
            fp.write("\n")
            line_number=1
            plt.xlim(0,12)
            plt.ylim(0,1)
            plt.xticks([])
            plt.gca().yaxis.set_major_formatter(ScalarFormatter())
            for key, value in infiles.items():
                plt.plot(line_number, res_avg_fraction_tot_integral[key][x,y,i], next(mmm), markersize=marker_size, color=next(colors), label=value[1])
                line_number=line_number+1
            plt.minorticks_on()
            plt.grid(b=True, color='#BBBBBB', linestyle='-')
            plt.legend(loc='lower right')
            plt.tight_layout()
            plt.savefig(od + "/4Dvolume_ratio_of_integrals_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".png", dpi=300, pad_inches=0.)
            plt.savefig(od + "/4Dvolume_ratio_of_integrals_system_comparison_Xlim_" + Xval + "_Ylim_" + Yval + "_edens_"+epsval+".pdf", pad_inches=0.)
            plt.close('all')
            fp.close()
