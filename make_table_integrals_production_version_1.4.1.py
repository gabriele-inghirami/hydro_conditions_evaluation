# version 1.4.1 - 20/12/2021 this program makes tables with integrals of the results

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import gzip

# minimum pressure density in GeV/fm^3
pmin = 0.0001

# we get the name output directory
N_input_files = len(sys.argv) - 1

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_table_integrals.py <pickled list of input files> <output directory>')
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

results = {}
errors = {}
results_avg = {}


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

    dt=tt[1]-tt[0]

    results[key] = N_cells_avg.sum(axis=0)*dt
    errors[key] = N_cells_sd.sum(axis=0)*dt

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    num=3.*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)

    n_eps = len(elim)
    n_X = len(Xlim)
    n_Y = len(Ylim)
    n_XT = len(XT_lim)
    n_RP = len(RP_lim)

    nx, ny, nz = lattice["dimensions"]

    counts=np.zeros((n_eps,n_X,n_Y),dtype=np.float64)
    for ie in range(n_eps):
        for ix in range(n_X):
            for iy in range(n_Y):
                for it in range(len(tt)):
                    for i in range(nx):
                        for j in range(ny):
                            for k in range(nz):
                                if(Tmunu[it,iT00,i,j,k]>elim[ie] and den[it,i,j,k]>=pmin and X_avg[it,i,j,k]<Xlim[ix] and Y_avg[it,i,j,k]<Ylim[iy]):
                                    counts[ie,ix,iy]=counts[ie,ix,iy]+1.
  
    counts=counts*dt
    results_avg[key] = counts.copy()
    counts = None
    data = None
    lattice = None
    tt = None
    N_events = None
    en_per_part = None
    T_munu = None
    X = None
    Y = None
    XT = None
    RP = None
    N_cells_avg = None
    N_cells_sd = None
    Non_empty_cells_avg = None
    Non_empty_cells_sd = None
    N_cells_avg_fraction = None
    N_cells_sd_fraction = None
    N_clumps_avg = None
    N_clumps_percent_avg = None
    N_clumps_single_cell_percent_tot = None
    N_clumps_single_cells_tot = None
    N_clumps_empty_tot = None


for x in range(n_X):
    Xval = '{:3.2f}'.format(Xlim[x])
    for y in range(n_Y):
                Yval = '{:3.2f}'.format(Ylim[y])
                fout=open(od + "/table_time_integrated_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                fout.write("# constraints: Xlim: " + Xval +", Ylim: " + Yval + "\n")
                fout.write("# first colum: reaction, remaining columns: volume in fm^4 above the energy density limit, with uncertainty\n")
                for i in range(n_eps):
                    fout.write("# column "+str(2+i*2)+": energy limit: "+'{:6.0f}'.format(elim[i] * 1000)+", column "+str(3+i*2)+": std deviation\n")
                fout.write("\n")     
                for key, value in results.items():
                    fout.write(key)
                    for i in range(n_eps):
                        fout.write("    "+'{:6.2f}'.format(value[i,x,y,0,0])+" ("+'{:6.2f}'.format(errors[key][i,x,y,0,0])+")")
                    fout.write("\n")     
                fout.close()

                fout=open(od + "/table_time_integrated_average_Tmunu_" + "_Xlim_" + Xval + "_Ylim_" + Yval + ".dat","w")
                fout.write("# constraints: Xlim: " + Xval +", Ylim: " + Yval + "\n")
                fout.write("# first colum: reaction, remaining columns: volume in fm^4 above the energy density limit\n")
                for i in range(n_eps):
                    fout.write("# column "+str(2+i*2)+"  energy limit: "+'{:6.0f}'.format(elim[i] * 1000)+" [MeV/fm^3]\n")
                fout.write("\n")     
                for key, value in results.items():
                    fout.write(key)
                    for i in range(n_eps):
                        fout.write("    "+'{:6.2f}'.format(results_avg[key][i,x,y]))
                    fout.write("\n")     
                fout.close()
