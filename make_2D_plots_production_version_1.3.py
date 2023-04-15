# version 1.3 - 10/12/2021 this program makes 2D plots of energy density vs X and Y 

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

matplotlib.use('Agg')

# bins for the histograms
#ebins = np.linspace(0.01,0.8,50,endpoint=True)
ebins = np.linspace(0.01,1.2,80,endpoint=True)

xbins = np.linspace(0.001,1.2,80,endpoint=True)
ybins = np.linspace(0.001,1.2,80,endpoint=True)

xbins_avg = np.linspace(0.001,1.2,80,endpoint=True)
ybins_avg = np.linspace(0.001,1.2,80,endpoint=True)

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

# we get the name output directory
N_input_files = len(sys.argv) - 1

# minimum value of T^11+T^22+T^33 to accept a grid point (units: GeV/fm^3)
pmin=0.0001

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_2D_plots.py <pickled list of input files> <output directory>')
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

# parameters for the plots
plt.rcParams.update({'font.size': 12})


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

    nt = len(tt)
    eps=Tmunu[:,iT00,:,:,:]
    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    num=3.*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)

    nx,ny,nz=lattice["dimensions"]

    #label_std=", event by event average"
    #label_avg=r", $T^{\mu\nu}$"+" average"

    # option to get normalized (True) or non normalized (False) histograms
    normalized=False

    # choice of the colormap
    #colmap = 'RdGy_r'
    colmap = 'inferno'

    # now we make the plots
    for i in range(nt):
        # we remove the points with energy density below the minimum threshold
        for ik in range(nz):
            for ij in range(ny):
                for ii in range(nx):
                    if (eps[i,ii,ij,ik]<elim[0]):
                        X[i,ii,ij,ik]=0.
                        Y[i,ii,ij,ik]=0.
                        X_avg[i,ii,ij,ik]=0.
                        Y_avg[i,ii,ij,ik]=0.
                        

        tstring = '{:3.2f}'.format(tt[i])+"fm"

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(X[i].flatten(), eps[i].flatten(), bins=(xbins, ebins), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_edens_vs_X_ebe_2Dhist_t_"+tstring+".dat","w")
        fp.write("# X_ebe      energy density [GeV/fm^3]        counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_std+", t="+tstring)
        img = plt.imshow(np.log10(H+0.1), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("X")
        plt.ylabel(r"$\varepsilon$ $[GeV/fm^3]$")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_edens_vs_X_ebe_2Dhist_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(X[i].flatten(), Y[i].flatten(), bins=(xbins, ybins), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_Y_ebe_vs_X_ebe_2Dhist_t_"+tstring+".dat","w")
        fp.write("# X      Y       counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_std+", t="+tstring)
        img = plt.imshow(np.log10(H+0.1), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_Y_ebe_vs_X_ebe_2Dhist_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(Y[i].flatten(), eps[i].flatten(), bins=(ybins, ebins), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_edens_vs_Y_ebe_2Dhist_t_"+tstring+".dat","w")
        fp.write("# Y      energy density [GeV/fm^3]        counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_std+", t="+tstring)
        img=plt.imshow(np.log10(H+0.01), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("Y")
        plt.ylabel(r"$\varepsilon$ $[GeV/fm^3]$")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_edens_vs_Y_ebe_2Dhist_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(X_avg[i].flatten(), eps[i].flatten(), bins=(xbins_avg, ebins), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_edens_vs_X_2Dhist_avg_t_"+tstring+".dat","w")
        fp.write("# X      energy density [GeV/fm^3]        counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_avg+", t="+tstring)
        img = plt.imshow(np.log10(H+0.01), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("X")
        plt.ylabel(r"$\varepsilon$ $[GeV/fm^3]$")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_edens_vs_X_2Dhist_avg_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(X_avg[i].flatten(), Y_avg[i].flatten(), bins=(xbins_avg, ybins_avg), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_Y_vs_X_2Dhist_avg_t_"+tstring+".dat","w")
        fp.write("# X      Y       counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_avg+", t="+tstring)
        img = plt.imshow(np.log10(H+0.01), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_Y_vs_X_2Dhist_avg_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')

        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 6
        fig_size[1] = 5
        plt.rcParams["figure.figsize"] = fig_size

        H, x_edges, y_edges = np.histogram2d(Y_avg[i].flatten(), eps[i].flatten(), bins=(ybins_avg, ebins), density=normalized)
        H = H.T
        fp=open(od + "/"+key+"_edens_vs_Y_2Dhist_avg_t_"+tstring+".dat","w")
        fp.write("# Y      energy density [GeV/fm^3]        counts\n")    
        xpoints = (x_edges[1:]+x_edges[:-1])/2.
        ypoints = (y_edges[1:]+y_edges[:-1])/2.
        for ix in range(len(xpoints)):
            for iy in range(len(ypoints)):
                fp.write('{:7.3f}'.format(xpoints[ix])+"      "+'{:7.3f}'.format(ypoints[iy])+"      "+'{:9.6e}'.format(H[iy,ix])+"\n")
            fp.write("\n")
        fp.close()
        #fig = plt.figure(figsize=(5, 3))
        #plt.title(value[1]+label_avg+", t="+tstring)
        img=plt.imshow(np.log10(H+0.01), interpolation='nearest', origin='lower', extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]], cmap=colmap)
        plt.clim(0,)
        plt.xlabel("Y")
        plt.ylabel(r"$\varepsilon$ $[GeV/fm^3]$")
        plt.colorbar(img,shrink=0.65,label="Log$_{10}$ counts")
        plt.tight_layout()
        plt.savefig(od + "/"+key+"_edens_vs_Y_2Dhist_avg_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')
