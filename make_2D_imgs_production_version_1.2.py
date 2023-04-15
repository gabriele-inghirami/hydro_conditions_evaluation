# version 1.2 - 10/12/2021 this program makes 2D plots of energy density, X and Y at different times
#                            for planes passing through the center of the grid orthogonal to one axis
#                            in the same figure epsilon, X, Y ebe, X, Y avgT



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

# minimum pressure for considering a certain cell part of the system (in GeV/fm^3)
pmin = 1.e-4

# we get the name output directory
N_input_files = len(sys.argv) - 1

# the pickled list of input files must be created with create_pickled_list_of_files_to_analize.py
if N_input_files != 2:
    print('Syntax: ./make_2D_imgs.py <pickled list of input files> <output directory>')
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

# colormap
chosen_cmp="hot"

# max and min values for the colormap
# average from event by event
hval_X=2
lval_X=0
hval_Y=2.0
lval_Y=0.
# average from average Tmunu
hval_XA=0.6
lval_XA=0
hval_YA=0.4
lval_YA=0.

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
    xp, yp, zp = lattice["dimensions"]
    dx, dy, dz = lattice["spacing"]
    ox, oy, oz = lattice["origin"]

    xx = np.linspace(ox, ox+xp*dx, num=xp, endpoint=True)
    yy = np.linspace(oy, oy+yp*dy, num=yp, endpoint=True)
    zz = np.linspace(oz, oz+zp*dz, num=zp, endpoint=True)

    # just for convenience and convention
    nx = xp
    ny = yp
    nz = zp

    i0=int(math.floor(xp/2))
    j0=int(math.floor(yp/2))
    k0=int(math.floor(zp/2))

    ps=25
    pe=-25

    den=Tmunu[:,iT11,:,:,:]+Tmunu[:,iT22,:,:,:]+Tmunu[:,iT33,:,:,:]
    num=abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT22,:,:,:])+abs(Tmunu[:,iT11,:,:,:]-Tmunu[:,iT33,:,:,:])+abs(Tmunu[:,iT33,:,:,:]-Tmunu[:,iT22,:,:,:])
    X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)
    num=3.*(abs(Tmunu[:,iT12,:,:,:])+abs(Tmunu[:,iT13,:,:,:])+abs(Tmunu[:,iT23,:,:,:]))
    Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den>=pmin)

    for it in range(nt):
        tstring='{:4.2f}'.format(tt[it])
        tnstring='{:04.2f}'.format(tt[it])
        #tstep='{:04d}'.format(it)

        fig_size = plt.rcParams["figure.figsize"]
        plt.rcParams.update({'font.size': 10})
        fig_size[0] = 12
        fig_size[1] = 4.
        plt.rcParams["figure.figsize"]=fig_size

        #plt.suptitle(value[1]+", t="+tstring+"fm",y=0.995,fontsize=6)

        fp=open(od + "/"+key+"_energy_density_t_"+tnstring+".dat","w")
        fp.write("# time: "+tstring+"\n")

        maxvalue=np.amax(Tmunu[it,iT00,ps:pe,ps:pe,k0])
        minvalue=np.amin(Tmunu[it,iT00,ps:pe,ps:pe,k0])
        plt.subplot(1,3,1)
        plt.imshow(Tmunu[it,iT00,ps:pe,ps:pe,k0].transpose(),extent=[xx[ps], xx[pe], yy[ps],yy[pe]], origin='lower',cmap=chosen_cmp,vmin=minvalue,vmax=maxvalue)
        plt.title(r'$\varepsilon$'+" [GeV/fm$^3$] at z=0")
        plt.xlabel('x [fm]')
        plt.ylabel('y [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# z = "+'{:6.3f}'.format(zz[k0])+"fm \n")
        fp.write("\n# x [fm]      y [fm]     en. dens. [GeV/fm^3]\n")
        for xq in range(ps,nx+pe):
            for yq in range(ps,ny+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(yy[yq])+"  "+'{:12.9e}'.format(Tmunu[it,iT00,xq,yq,k0])+"\n")
            fp.write("\n")
        fp.write("\n")
                    
        maxvalue=np.amax(Tmunu[it,iT00,ps:pe,j0,ps:pe])
        minvalue=np.amin(Tmunu[it,iT00,ps:pe,j0,ps:pe])
        plt.subplot(1,3,2)
        plt.imshow(Tmunu[it,iT00,ps:pe,j0,ps:pe].transpose(),extent=[xx[ps], xx[pe], zz[ps],zz[pe]], origin='lower',cmap=chosen_cmp,vmin=minvalue,vmax=maxvalue)
        plt.title(r'$\varepsilon$'+" [GeV/fm$^3$] at y=0")
        plt.xlabel('x [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# y = "+'{:6.3f}'.format(yy[j0])+"fm \n")
        fp.write("\n# x [fm]      z [fm]     en. dens.T^00 [GeV/fm^3]\n")
        for xq in range(ps,nx+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Tmunu[it,iT00,xq,j0,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        maxvalue=np.amax(Tmunu[it,iT00,i0,ps:pe,ps:pe])
        minvalue=np.amin(Tmunu[it,iT00,i0,ps:pe,ps:pe])
        plt.subplot(1,3,3)
        plt.imshow(Tmunu[it,iT00,i0,ps:pe,ps:pe].transpose(),extent=[yy[ps], yy[pe], zz[ps],zz[pe]], origin='lower',cmap=chosen_cmp,vmin=minvalue,vmax=maxvalue)
        plt.title(r'$\varepsilon$'+" [GeV/fm$^3$] at x=0")
        plt.xlabel('y [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# x = "+'{:6.3f}'.format(xx[i0])+"fm \n")
        fp.write("\n# y [fm]      z [fm]     en. dens.T^00 [GeV/fm^3]\n")
        for yq in range(ps,ny+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(yy[yq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Tmunu[it,iT00,i0,yq,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        fp.close()

        plt.tight_layout()

        plt.savefig(od + "/"+key+"_energy_density_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')




        fp=open(od + "/"+key+"_X_t_"+tnstring+".dat","w")
        fp.write("# time: "+tstring+"\n")

        fig_size = plt.rcParams["figure.figsize"]
        plt.rcParams.update({'font.size': 10})
        fig_size[0] = 12
        fig_size[1] = 8
        plt.rcParams["figure.figsize"]=fig_size

        #plt.suptitle(value[1]+", t="+tstring+"fm",y=0.995,fontsize=6)

        plt.subplot(2,3,1)
        plt.imshow(X[it,ps:pe,ps:pe,k0].transpose(),extent=[xx[ps], xx[pe], yy[ps], yy[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_X,vmax=hval_X)
        plt.title("$X_{ebe}$ at z=0")
        plt.xlabel('x [fm]')
        plt.ylabel('y [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# z = "+'{:6.3f}'.format(zz[k0])+"fm \n")
        fp.write("\n# x [fm]      y [fm]     X_ebe (average value of pressure anisotropy from many single events)\n")
        for xq in range(ps,nx+pe):
            for yq in range(ps,ny+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(yy[yq])+"  "+'{:12.9e}'.format(X[it,xq,yq,k0])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,2)
        plt.imshow(X[it,ps:pe,j0,ps:pe].transpose(),extent=[xx[ps], xx[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_X,vmax=hval_X)
        plt.title("$X_{ebe}$ at y=0")
        plt.xlabel('x [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# y = "+'{:6.3f}'.format(yy[j0])+"fm \n")
        fp.write("\n# x [fm]      z [fm]     X_ebe (average value of pressure anisotropy from many single events)\n")
        for xq in range(ps,nx+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(X[it,xq,j0,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,3)
        plt.imshow(X[it,i0,ps:pe,ps:pe].transpose(),extent=[yy[ps], yy[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_X,vmax=hval_X)
        plt.title("$X_{ebe}$ at x=0")
        plt.xlabel('y [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# x = "+'{:6.3f}'.format(xx[i0])+"fm \n")
        fp.write("\n# y [fm]      z [fm]     X_ebe (average value of pressure anisotropy from many single events)\n")
        for yq in range(ps,ny+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(yy[yq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(X[it,i0,yq,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,4)
        plt.imshow(X_avg[it,ps:pe,ps:pe,k0].transpose(),extent=[xx[ps], xx[pe], yy[ps], yy[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_XA,vmax=hval_XA)
        plt.title("$X$ at z=0")
        plt.xlabel('x [fm]')
        plt.ylabel('y [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# z = "+'{:6.3f}'.format(zz[k0])+"fm \n")
        fp.write("\n# x [fm]      y [fm]     X (pressure anisotropy of average energy momentum tensor)\n")
        for xq in range(ps,nx+pe):
            for yq in range(ps,ny+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(yy[yq])+"  "+'{:12.9e}'.format(X[it,xq,yq,k0])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,5)
        plt.imshow(X_avg[it,ps:pe,j0,ps:pe].transpose(),extent=[xx[ps], xx[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_XA,vmax=hval_XA)
        plt.title("$X$ at y=0")
        plt.xlabel('x [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# y = "+'{:6.3f}'.format(yy[j0])+"fm \n")
        fp.write("\n# x [fm]      z [fm]     X (pressure anisotropy of average energy momentum tensor)\n")
        for xq in range(ps,nx+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(X[it,xq,j0,zq])+"\n")
            fp.write("\n")
        fp.write("\n")
        
        plt.subplot(2,3,6)
        plt.imshow(X_avg[it,i0,ps:pe,ps:pe].transpose(),extent=[yy[ps], yy[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_XA,vmax=hval_XA)
        plt.title("$X$ at x=0")
        plt.xlabel('y [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# x = "+'{:6.3f}'.format(xx[i0])+"fm \n")
        fp.write("\n# y [fm]      z [fm]     X (pressure anisotropy of average energy momentum tensor)\n")
        for yq in range(ps,ny+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(yy[yq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(X[it,i0,yq,zq])+"\n")
            fp.write("\n")
        fp.write("\n")
        
        fp.close()

        plt.tight_layout()

        plt.savefig(od + "/"+key+"_X_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')



        fp=open(od + "/"+key+"_Y_t_"+tnstring+".dat","w")
        fp.write("# time: "+tstring+"\n")


        fig_size = plt.rcParams["figure.figsize"]
        plt.rcParams.update({'font.size': 10})
        fig_size[0] = 12
        fig_size[1] = 8
        plt.rcParams["figure.figsize"]=fig_size

        #plt.suptitle(value[1]+", t="+tstring+"fm",y=0.995,fontsize=6)

        plt.subplot(2,3,1)
        plt.imshow(Y[it,ps:pe,ps:pe,k0].transpose(),extent=[xx[ps], xx[pe], yy[ps], yy[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_Y,vmax=hval_Y)
        plt.title("$Y_{ebe}$ at z=0")
        plt.xlabel('x [fm]')
        plt.ylabel('y [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# z = "+'{:6.3f}'.format(zz[k0])+"fm \n")
        fp.write("\n# x [fm]      y [fm]     Y_ebe (average value of off-diagonality from many single events)\n")
        for xq in range(ps,nx+pe):
            for yq in range(ps,ny+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(yy[yq])+"  "+'{:12.9e}'.format(Y[it,xq,yq,k0])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,2)
        plt.imshow(Y[it,ps:pe,j0,ps:pe].transpose(),extent=[xx[ps], xx[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_Y,vmax=hval_Y)
        plt.title("$Y_{ebe}$ at y=0")
        plt.xlabel('x [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# y = "+'{:6.3f}'.format(yy[j0])+"fm \n")
        fp.write("\n# x [fm]      z [fm]     Y_ebe (average value of off-diagonality from many single events)\n")
        for xq in range(ps,nx+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Y[it,xq,j0,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,3)
        plt.imshow(Y[it,i0,ps:pe,ps:pe].transpose(),extent=[yy[ps], yy[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_Y,vmax=hval_Y)
        plt.title("$Y_{ebe}$ at x=0")
        plt.xlabel('y [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# x = "+'{:6.3f}'.format(xx[i0])+"fm \n")
        fp.write("\n# y [fm]      z [fm]     Y_ebe (average value of off-diagonality from many single events)\n")
        for yq in range(ps,ny+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(yy[yq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Y[it,i0,yq,zq])+"\n")
            fp.write("\n")
        fp.write("\n")


        plt.subplot(2,3,4)
        plt.imshow(Y_avg[it,ps:pe,ps:pe,k0].transpose(),extent=[xx[ps], xx[pe], yy[ps], yy[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_YA,vmax=hval_YA)
        plt.title("$Y$ at z=0")
        plt.xlabel('x [fm]')
        plt.ylabel('y [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# z = "+'{:6.3f}'.format(zz[k0])+"fm \n")
        fp.write("\n# x [fm]      y [fm]     Y (off-diagonality of average energy momentum tensor)\n")
        for xq in range(ps,nx+pe):
            for yq in range(ps,ny+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(yy[yq])+"  "+'{:12.9e}'.format(Y[it,xq,yq,k0])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,5)
        plt.imshow(Y_avg[it,ps:pe,j0,ps:pe].transpose(),extent=[xx[ps], xx[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_YA,vmax=hval_YA)
        plt.title("$Y$ at y=0")
        plt.xlabel('x [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# y = "+'{:6.3f}'.format(yy[j0])+"fm \n")
        fp.write("\n# x [fm]      z [fm]     Y (off-diagonality of average energy momentum tensor)\n")
        for xq in range(ps,nx+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(xx[xq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Y[it,xq,j0,zq])+"\n")
            fp.write("\n")
        fp.write("\n")

        plt.subplot(2,3,6)
        plt.imshow(Y_avg[it,i0,ps:pe,ps:pe].transpose(),extent=[yy[ps], yy[pe], zz[ps], zz[pe]], origin='lower',cmap=chosen_cmp,vmin=lval_YA,vmax=hval_YA)
        plt.title("$Y$ at x=0")
        plt.xlabel('y [fm]')
        plt.ylabel('z [fm]')
        ax=plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax,format="%5.3f")
        fp.write("\n# x = "+'{:6.3f}'.format(xx[i0])+"fm \n")
        fp.write("\n# y [fm]      z [fm]     Y (off-diagonality of average energy momentum tensor)\n")
        for yq in range(ps,ny+pe):
            for zq in range(ps,nz+pe):
                fp.write('{:6.3f}'.format(yy[yq])+"  "+'{:6.3f}'.format(zz[zq])+"  "+'{:12.9e}'.format(Y[it,i0,yq,zq])+"\n")
            fp.write("\n")
        fp.write("\n")
        
        fp.close()


        #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.tight_layout()
        #plt.subplots_adjust(left=0.05,
        #            bottom=0.02, 
        #            right=0.95, 
        #            top=0.97, 
        #            wspace=0.25, 
        #            hspace=0.3)

        plt.savefig(od + "/"+key+"_Y_t_"+tstring+".png", dpi=300, pad_inches=0.)
        plt.close('all')
