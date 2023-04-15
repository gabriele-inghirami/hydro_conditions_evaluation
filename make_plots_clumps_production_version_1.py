# version 1.1.1 - 09/02/2022 this program plots the number of clumps vs time 

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

# we use a dictionary for the files to be analized
# the key will enter in the output filenames
# the values are a tuple with: 0) the name of the files 1) a string with the common part of the plot titles
infiles = {}
infiles["Au_Elab_1_23_central_standard"] = ("results_Au_Elab_1_23_wcl.pickle.gz", "Au, Elab=1.23AGeV")

# path of the directory containing the data
datapath = "./pickled_data_final/"

# colors of the points in the plots
colors=["orange","red","magenta"]

# end of the plot time range
time_end=30

# we check if the directory with the data exists
if not os.path.exists(datapath):
    print("Error, the directory " + datapath + " , that should contain the input data does not exist!")
    sys.exit(1)

# we get the name output directory
N_input_files = len(sys.argv) - 1

if N_input_files != 1:
    print('Syntax: ./make_plot_clumps.py <output directory>')
    sys.exit(1)

od = sys.argv[1]
if (not os.path.exists(od)):
    os.mkdir(od)

results = {}
errors = {}

symbs = ("o","v","^")

# actually we work with just item and the code implicitly assumes that we have only one
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
    N_cells_sd_fraction, N_clumps_avg,N_clumps_percent_tot,N_clumps_single_cell_percent_tot,N_clumps_single_cell_tot,N_clumps_nonempty_tot = data[:]

N_nonsingle_clumps=N_clumps_avg*N_events-N_clumps_single_cell_tot

n_eps = len(elim)
n_X = len(Xlim)
n_Y = len(Ylim)

tt=np.array(tt)
iend=np.argmin(abs(tt-time_end))

for x in range(n_X):
    Xval = '{:3.2f}'.format(Xlim[x])
    for y in range(n_Y):
                Yval = '{:3.2f}'.format(Ylim[y])
                fo=open(od + "/clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                fo.write("# Au+Au, E$_{lab}$=1.23 AGeV, b=0-3.3fm"+ ", X<" + Xval + ", Y<" + Yval + "\n")
                plt.xlabel('t [fm]')
                plt.ylabel('<number of clumps>')
                plt.xlim(-3, tt[iend])
                for ie in range(n_eps):
                    epsval = '{:3.0f}'.format(1000.*elim[ie])
                    fo.write("\n# Energy density > "+epsval+"\n")
                    fo.write("# time [fm]      <number of clumps>\n")
                    plt.plot(tt[:iend], N_clumps_avg[:iend,ie,x,y,0,0], symbs[ie], color=colors[ie], markersize=2, label=epsval+"MeV/fm$^3$")
                    for it in range(iend):
                        fo.write('{:5.2f}'.format(tt[it])+"      "+'{:9.5f}'.format(N_clumps_avg[it,ie,x,y,0,0])+"\n")
                    fo.write("\n")
                plt.minorticks_on()
                plt.grid(b=True, color='#BBBBBB', linestyle='-')
                plt.legend()
                plt.tight_layout()
                plt.savefig(od + "/clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".png", dpi=300, pad_inches=0.)
                plt.savefig(od + "/clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".pdf", pad_inches=0.)
                plt.close('all')
                fo.close()


                fo=open(od + "/compactness_of_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                fo.write("# Au+Au, E$_{lab}$=1.23 AGeV, b=0-3.3fm"+ ", X<" + Xval + ", Y<" + Yval + "\n")
                plt.xlabel('t [fm]')
                plt.ylabel('<compactness of clumps>')
                plt.xlim(-3, tt[iend])
                for ie in range(n_eps):
                    epsval = '{:3.0f}'.format(1000.*elim[ie])
                    fo.write("\n# Energy density > "+epsval+"\n")
                    fo.write("# time [fm]      <compactness>\n")
                    ratio_data=np.divide(N_clumps_percent_tot[:iend,ie,x,y,0,0],N_nonsingle_clumps[:iend,ie,x,y,0,0],out=np.zeros(N_clumps_percent_tot[:iend,ie,x,y,0,0].shape, dtype=np.float64), where=N_nonsingle_clumps[:iend,ie,x,y,0,0]!=0)
                    plt.plot(tt[:iend], ratio_data, symbs[ie], color=colors[ie], markersize=2, label=epsval+"MeV/fm$^3$")
                    for it in range(iend):
                        fo.write('{:5.2f}'.format(tt[it])+"      "+'{:7.5f}'.format(ratio_data[it])+"\n")
                    fo.write("\n")
                plt.minorticks_on()
                plt.grid(b=True, color='#BBBBBB', linestyle='-')
                plt.legend()
                plt.tight_layout()
                plt.savefig(od + "/compactness_of_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".pdf", pad_inches=0.)
                plt.savefig(od + "/compactness_of_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".png", dpi=300, pad_inches=0.)
                plt.close('all')
                fo.close()

                fo=open(od + "/fraction_of_single_cell_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                fo.write("# Au+Au, E$_{lab}$=1.23 AGeV, b=0-3.3fm"+ ", X<" + Xval + ", Y<" + Yval + "\n")
                plt.xlabel('t [fm]')
                plt.ylabel('<fraction of single cell clumps>')
                plt.xlim(-3, tt[iend])
                for ie in range(n_eps):
                    epsval = '{:3.0f}'.format(1000.*elim[ie])
                    fo.write("\n# Energy density > "+epsval+"\n")
                    fo.write("# time [fm]      <fraction of single cell clumps>\n")
                    ratio_data=np.divide(N_clumps_single_cell_percent_tot[:iend,ie,x,y,0,0],N_clumps_nonempty_tot[:iend,ie,x,y,0,0],out=np.zeros(N_clumps_single_cell_percent_tot[:iend,ie,x,y,0,0].shape, dtype=np.float64), where=N_clumps_nonempty_tot[:iend,ie,x,y,0,0]!=0)
                    plt.plot(tt[:iend], ratio_data, symbs[ie], color=colors[ie], markersize=2, label=epsval+"MeV/fm$^3$")
                    for it in range(iend):
                        fo.write('{:5.2f}'.format(tt[it])+"      "+'{:7.5f}'.format(ratio_data[it])+"\n")
                    fo.write("\n")
                plt.minorticks_on()
                plt.grid(b=True, color='#BBBBBB', linestyle='-')
                plt.legend()
                plt.tight_layout()
                plt.savefig(od + "/fraction_of_single_cell_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".png", dpi=300, pad_inches=0.)
                plt.savefig(od + "/fraction_of_single_cell_clumps_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".pdf", pad_inches=0.)
                plt.close('all')
                fo.close()
