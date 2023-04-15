# analize_data.py - version 0.8.5 - 15/10/2021

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import os.path
import glob
from timeit import default_timer as timer

sys.setrecursionlimit(10000)

#use binary (True) or ascii (False) input data
use_binary=True

#density tpye (it can be baryon or hadron)
density_type="hadron"
#density_type="baryon"

#if False it prints only error messages, if True it writes what it is doing at the moment and the intermediate results 
verbose=False

#if True it searches for clumps
search_clumps=True

#output file version check
designed_lattice_version=1.0

#if True it prints only the energy density time evolution in the point with given lattice indexes
print_only_energy_density_evolution=False

#pressure threshold to consider a cell non empty (in GeV/fm^3)
p_min_for_non_emptyness=0.0001

#energy density thresholds (in GeV/fm^3) to count cells (tuple)
#elim=(0.1,0.2,0.3,0.4,0.5)
elim=(0.001,0.1,0.5)

#pressure anisotropy thresholds
#Xlim=(1.0,0.75,0.5,0.25,0.1)
Xlim=(0.5,0.3,0.2,0.1)

#off-diagonality thresholds
#Ylim=(1.0,0.75,0.5,0.25,0.1)
Ylim=(0.5,0.3,0.2,0.1)

#transverse pressure anisotropy thresholds
#XT_lim=(0.75,0.5,0.25,0.1)
XT_lim=(10000.,)

#ratio between transverse and longitudinal pressure
#RP_lim=(1.,0.5,0.25,0.1)
RP_lim=(10000.,)

#timing format
tf='{:8.4f}'

if(verbose):
    init_start=timer()

#we parse the command line arguments
N_input_args=len(sys.argv)-1

if(N_input_args!=2):
   print ('Syntax: ./analize_data.py <data dir> <outputfile>')
   print ("Directory containing the files produced by SMASH with Thermodynamic Lattice Output")
   print ("outputfile is obviously the name of the output file with the results of the postprocessing")
   sys.exit(1)

#we get the name of input and output files
inputdir=sys.argv[1]
outputfile=sys.argv[2]

#we prepare lists of the input files
if (density_type == "hadron"):
    if(use_binary):
        net_bar_files=glob.glob(inputdir+'/hadron_j_QBS_*.bin')
        Tmunu_files=glob.glob(inputdir+'/hadron_tmn_landau_*.bin')
        vLandau_files=glob.glob(inputdir+'/hadron_v_landau_*.bin')
    else:
        net_bar_files=glob.glob(inputdir+'/hadron_j_QBS_*.dat')
        Tmunu_files=glob.glob(inputdir+'/hadron_tmn_landau_*.dat')
        vLandau_files=glob.glob(inputdir+'/hadron_v_landau_*.dat')
elif (density_type == "baryon"):
    if(use_binary):
        net_bar_files=glob.glob(inputdir+'/net_baryon_j_QBS_*.bin')
        Tmunu_files=glob.glob(inputdir+'/net_baryon_tmn_landau_*.bin')
        vLandau_files=glob.glob(inputdir+'/net_baryon_v_landau_*.bin')
    else:
        net_bar_files=glob.glob(inputdir+'/net_baryon_j_QBS_*.dat')
        Tmunu_files=glob.glob(inputdir+'/net_baryon_tmn_landau_*.dat')
        vLandau_files=glob.glob(inputdir+'/net_baryon_v_landau_*.dat')
else:
    print("Unknown density_type parameter (please, check the first lines of the script source code and fix it)")
    sys.exit(2)

#we sort the lists of input files
net_bar_files.sort()
Tmunu_files.sort()
vLandau_files.sort()

nf_bar=len(net_bar_files)
nf_tmn=len(Tmunu_files)
nf_vl=len(vLandau_files)

if((nf_bar != nf_tmn) or (nf_bar != nf_vl)):
    print("Sorry, but I can't continue.")
    print("I have found "+str(nf_bar)+" density current files, "+str(nf_tmn)+" Tmunu files, "+str(nf_vl)+" Landau velocity files")
    sys.exit(2)

#we check that all the input files have non zero length
for i in range(nf_bar):
    if((os.path.getsize(net_bar_files[i])==0) or (os.path.getsize(Tmunu_files[i])==0) or (os.path.getsize(vLandau_files[i])==0)):
        print("Because of a zero length file, I will not consider:")
        print(nf_bar.pop(i))
        print(nf_tmn.pop(i))
        print(nf_vl.pop(i))

#we update the length of the files
nf_bar=len(net_bar_files)
nf_tmn=len(Tmunu_files)
nf_vl=len(vLandau_files)

if(nf_bar*nf_tmn*nf_vl==0):
    print("Input files missing. I quit")
    sys.exit(2)
else:
    nf=nf_bar #we use a common variable for all file lengths

#dictionary containing information about the grid
lattice={}

#number of events
N_events=0

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

#indexes corresponding to jB[0:3] in j_QBS
kb0=4
kb1=5
kb2=6
kb3=7

#empty list with the results
T=[]
en_per_part=[]
X=[]
Y=[]
XT=[]
RP=[]
N_clumps=[]
N_clumps_single_cell=[]
N_clumps_percent=[]
N_clumps_single_cell_percent=[]
N_clumps_nonempty=[]
N_cells_series=[]
N_cells_series_fraction=[]
Non_empty_cells_series=[]

#empty list with the output times
tt=[]

#info about cells in clumps
cl_cells=[]
cl_sizes=[]

#functions to read the header
def read_ascii_header(infile):
    lattice_version=float(infile.readline().split()[4])
    if(lattice_version != designed_lattice_version):
        print("Sorry, this code for analysis is designed for output version: "+str(designed_lattice_version)+", while you provided "+str(lattice_version)+"\n. I quit.")
        sys.exit(2)
    lattice_quantity=infile.readline().split()[1]
    lattice_dimensions=np.int32(infile.readline().split()[2:5])
    lattice_spacing=np.float64(infile.readline().split()[2:5])
    lattice_origin=np.float64(infile.readline().split()[2:5])
    return lattice_dimensions, lattice_spacing, lattice_origin

def read_binary_header(infile):
    lattice_version=np.fromfile(infile,dtype=np.float64,count=1)[0]
    if(lattice_version != designed_lattice_version):
        print("Sorry, this code for analysis is designed for output version: "+str(designed_lattice_version)+", while you provided "+str(lattice_version)+"\n. I quit.")
        sys.exit(2)
    lattice_quantity=np.fromfile(infile,dtype=np.int32,count=1)[0]
    lattice_dimensions=np.fromfile(infile,dtype=np.int32,count=3)
    lattice_spacing=np.fromfile(infile,dtype=np.float64,count=3)
    lattice_origin=np.fromfile(infile,dtype=np.float64,count=3)
    return lattice_dimensions, lattice_spacing, lattice_origin

def read_ascii_timestep_tmn(infile):
    #it returns True in case of EoF
    time_entry=infile.readline()
    if(time_entry==''):
        return True, None, None
    else:
        time=np.float64(time_entry)
    Tmunu=np.empty((10,nx,ny,nz),dtype=np.float64)
    for i in range(10):
        tmp_arr=np.loadtxt(infile,max_rows=ny*nz).flatten().reshape(nz,ny,nx).transpose()
        Tmunu[i,:,:,:]=tmp_arr.copy()
    return False, time, Tmunu

def read_binary_timestep_tmn(infile):
    #it returns True in case of EoF
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    if(len(time_entry)==0):
        return True, None, None
    else:
        time=time_entry[0]
    Tmunu=np.empty((10,nx,ny,nz),dtype=np.float64)
    for i in range(10):
        tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz).flatten().reshape(nz,ny,nx).transpose()
        Tmunu[i,:,:,:]=tmp_arr.copy()
    return False, time, Tmunu

def read_ascii_timestep_jqbs(infile):
    #we do not check the time entry
    time_entry=infile.readline()
    tmp_arr=np.loadtxt(infile,max_rows=nx*ny*nz).flatten().reshape(nz,ny,nx,12).transpose()
    return tmp_arr

def read_binary_timestep_jqbs(infile):
    #we do not check the time entry
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz*12).flatten().reshape(nz,ny,nx,12).transpose()
    return tmp_arr

def read_ascii_timestep_vl(infile):
    #we do not check the time entry
    time_entry=infile.readline()
    tmp_arr=np.loadtxt(infile,max_rows=nx*ny*nz).flatten().reshape(nz,ny,nx,3).transpose()
    return tmp_arr

def read_binary_timestep_vl(infile):
    #we do not check the time entry
    time_entry=np.fromfile(infile,dtype=np.float64,count=1)
    tmp_arr=np.fromfile(infile,dtype=np.float64,count=nx*ny*nz*3).flatten().reshape(nz,ny,nx,3).transpose()
    return tmp_arr

def find_clumps(ix,iy,iz,stack_depth):
    global ccc, cl_cells, cl_sizes
    if(verbose):
        print("Entering into: "+str(ix)+"  "+str(iy)+"  "+str(iz)+"  "+str(stack_depth))
    if(ccc[ix,iy,iz]!=0):
        if(verbose):
            print(str(ix)+"  "+str(iy)+"  "+str(iz)+" already checked")
        return
    if((denominator[ix,iy,iz]>=p_min_for_non_emptyness) and (Tmunu[iT00,ix,iy,iz]>0) and (Tmunu[iT00,ix,iy,iz]>elim_val) and (X_ratio[ix,iy,iz]<xlim_val) and (Y_ratio[ix,iy,iz]<ylim_val) and (XT_ratio[ix,iy,iz]<xt_lim_val) and (RP_entry[ix,iy,iz]<rp_lim_val)):
        if(verbose):
            print("Got: "+str(ix)+"  "+str(iy)+"  "+str(iz)+" SD: "+str(stack_depth))
        if(stack_depth==0): #if 0 it means that we are searching for a new clump
            cl_cells.append(1)
            cl_sizes.append(((ix,ix),(iy,iy),(iz,iz)))
            if(verbose):
                print("NEW")
                print(str(cl_cells))
        else:
            if(verbose):
                print("OLD")
                print(str(cl_cells))
            cl_cells[-1]=cl_cells[-1]+1
            cl_sizes[-1]=((min(ix,cl_sizes[-1][0][0]),max(ix,cl_sizes[-1][0][1])),(min(iy,cl_sizes[-1][1][0]),max(iy,cl_sizes[-1][1][1])),(min(iz,cl_sizes[-1][2][0]),max(iz,cl_sizes[-1][2][1])))
        ccc[ix,iy,iz]=len(cl_cells)
        for i in range(max(ix-1,0),min(ix+2,nx)):
         for j in range(max(iy-1,0),min(iy+2,ny)):
          for k in range(max(iz-1,0),min(iz+2,nz)):
              if(verbose):
                  print("Going deeper between: "+str(max(ix-1,0))+"  "+str(min(ix+1,nx))+"  "+str(max(iy-1,0))+"  "+str(min(iy+1,ny))+"  "+str(max(iz-1,0))+"  "+str(min(iz+1,nz)))
              find_clumps(i,j,k,stack_depth+1)
    else:
        if(verbose):
            print(str(ix)+"  "+str(iy)+"  "+str(iz)+" not included")
        ccc[ix,iy,iz]=-1
          

for n_i in range(nf):
    i_tmn=Tmunu_files[n_i]
    i_jqbs=net_bar_files[n_i]
    i_vl=vLandau_files[n_i]
    if(verbose):
        print("Opening "+i_tmn)
        start_time = timer()
    if(use_binary):
        if(i_tmn[-3:]==".gz"):
            fp_tmn=gzip.open(i_tmn,"rb")
        else:
            fp_tmn=open(i_tmn,"rb")
        lattice_dimensions, lattice_spacing, lattice_origin = read_binary_header(fp_tmn)
    else:
        if(i_tmn[-3:]==".gz"):
            fp_tmn=gzip.open(i_tmn,"r")
        else:
            fp_tmn=open(i_tmn,"r")
        lattice_dimensions, lattice_spacing, lattice_origin = read_ascii_header(fp_tmn)
    if(verbose):
        print("Opening "+i_jqbs)
        start_time = timer()
    if(use_binary):
        if(i_jqbs[-3:]==".gz"):
            fp_jqbs=gzip.open(i_jqbs,"rb")
        else:
            fp_jqbs=open(i_jqbs,"rb")
        lattice_dimensions_jqbs, lattice_spacing_jqbs, lattice_origin_jqbs = read_binary_header(fp_jqbs)
    else:
        if(i_jqbs[-3:]==".gz"):
            fp_jqbs=gzip.open(i_jqbs,"r")
        else:
            fp_jqbs=open(i_jqbs,"r")
        lattice_dimensions_jqbs, lattice_spacing_jqbs, lattice_origin_jqbs = read_ascii_header(fp_jqbs)
    if(verbose):
        print("Opening "+i_vl)
        start_time = timer()
    if(use_binary):
        if(i_vl[-3:]==".gz"):
            fp_vl=gzip.open(i_vl,"rb")
        else:
            fp_vl=open(i_vl,"rb")
        lattice_dimensions_vl, lattice_spacing_vl, lattice_origin_vl = read_binary_header(fp_vl)
    else:
        if(i_vl[-3:]==".gz"):
            fp_vl=gzip.open(i_vl,"r")
        else:
            fp_vl=open(i_vl,"r")
        lattice_dimensions_vl, lattice_spacing_vl, lattice_origin_vl = read_ascii_header(fp_vl)

    #we skip the check that the grid data for the various quantities are compatible with each other 

    if(n_i==0): #initially we need to acquire some information
        lattice["dimensions"]=lattice_dimensions
        lattice["spacing"]=lattice_spacing
        lattice["origin"]=lattice_origin
        nx,ny,nz=lattice_dimensions[:]
    else:
        if(not np.array_equal(lattice["dimensions"],lattice_dimensions)):
            print("Error in file "+i_tmn+": different lattice dimensions. Until now: "+str(lattice["dimensions"])+", this time: "+str(lattice_dimensions)+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice["spacing"],lattice_spacing)):
            print("Error in file "+i_tmn+": different lattice spacing. Until now: "+str(lattice["spacing"])+", this time: "+str(lattice_spacing)+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice["origin"],lattice_origin)):
            print("Error in file "+i_tmn+": different lattice origin. Until now: "+str(lattice["origin"])+", this time: "+str(lattice_origin)+".\nI quit.")
            sys.exit(2)

    index=0
    N_cells_event=[]
    N_cells_event_fraction=[]
    Non_empty_cells_event=[]
    N_clumps_event=[]
    N_clumps_single_cell_event=[]
    N_clumps_percent_event=[]
    N_clumps_single_cell_percent_event=[]
    N_clumps_nonempty_event=[]

    while(True):
        if(use_binary):
            eof,time,Tmunu=read_binary_timestep_tmn(fp_tmn)
            if(eof):
                break
            j_QBS=read_binary_timestep_jqbs(fp_jqbs) #first index: j component then x, y, z
            vl=read_binary_timestep_vl(fp_vl) #first index: v component then x, y, z
        else:
            eof,time,Tmunu=read_ascii_timestep_tmn(fp_tmn)
            if(eof):
                break
            j_QBS=read_ascii_timestep_jqbs(fp_jqbs) #first index: j component then x, y, z
            vl=read_ascii_timestep_vl(fp_vl) #first index: v component then x, y, z
        if(n_i==0):
            tt.append(time)
        else:
            if(time!=tt[index]):
                print("Error when reading file "+i+":")
                print("At step "+str(index)+" time "+str(tt[index])+" was expected, while "+str(time)+" was found. I quit.\n")
                sys.exit(2)
        glf=1./np.sqrt(1-vl[0,:,:,:]**2-vl[1,:,:,:]**2-vl[2,:,:,:]**2)
        ul_cov=(glf,-glf*vl[0,:,:,:],-glf*vl[1,:,:,:],-glf*vl[2,:,:,:])
        # the entries in the first dimension of j_QBS are: electric charge current, baryon current, strangeness current
        rho_B_Landau=j_QBS[4,:,:,:]*ul_cov[0]+j_QBS[5,:,:,:]*ul_cov[1]+j_QBS[6,:,:,:]*ul_cov[2]+j_QBS[7,:,:,:]*ul_cov[3]
        denominator=Tmunu[iT11,:,:,:]+Tmunu[iT22,:,:,:]+Tmunu[iT33,:,:,:]
        X_numerator=np.abs(Tmunu[iT11,:,:,:]-Tmunu[iT22,:,:,:])+np.abs(Tmunu[iT22,:,:,:]-Tmunu[iT33,:,:,:])+np.abs(Tmunu[iT33,:,:,:]-Tmunu[iT11,:,:,:])
        Y_numerator=3*(np.abs(Tmunu[iT12,:,:,:])+np.abs(Tmunu[iT23,:,:,:])+np.abs(Tmunu[iT13,:,:,:]))
        energy_per_hadron=np.divide(Tmunu[iT00,:,:,:],rho_B_Landau,out=np.zeros(Tmunu[iT00,:,:,:].shape, dtype=np.float64), where=rho_B_Landau!=0)
        X_ratio=np.divide(X_numerator,denominator,out=np.zeros(X_numerator.shape, dtype=np.float64), where=denominator!=0)
        Y_ratio=np.divide(Y_numerator,denominator,out=np.zeros(Y_numerator.shape, dtype=np.float64), where=denominator!=0)
        numerator_xy=np.abs(Tmunu[iT11,:,:,:]-Tmunu[iT22,:,:,:])
        sum_xy=Tmunu[iT11,:,:,:]+Tmunu[iT22,:,:,:]
        XT_ratio=np.divide(numerator_xy,sum_xy,out=np.zeros(numerator_xy.shape, dtype=np.float64), where=sum_xy!=0)
        PT_PR_ratio=np.divide(sum_xy,Tmunu[iT33,:,:,:],out=np.zeros(sum_xy.shape, dtype=np.float64), where=Tmunu[iT33,:,:,:]!=0)
        RP_entry=np.abs(0.5*PT_PR_ratio-1) 
        N_cells_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.int32)
        Non_empty_cells=0

        if(search_clumps):
            N_clumps_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.int32)
            N_clumps_single_cell_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.int32)
            N_clumps_nonempty_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.int32)
            N_clumps_percent_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
            N_clumps_single_cell_percent_entry=np.zeros((len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
            for ie in range(len(elim)):
             for xl in range(len(Xlim)):
              for yl in range(len(Ylim)):
               for tl in range(len(XT_lim)):
                for pl in range(len(RP_lim)):
                    ccc=np.zeros((nx,ny,nz),dtype=np.int32)
                    cl_cells=[]
                    cl_sizes=[]

                    elim_val=elim[ie]
                    xlim_val=Xlim[xl]
                    ylim_val=Ylim[yl]
                    xt_lim_val=XT_lim[tl]
                    rp_lim_val=RP_lim[pl]

                    if(verbose):
                        print("Checking clumps for: "+str(elim_val)+"  "+str(xlim_val)+"  "+str(ylim_val)+"  "+str(xt_lim_val)+"  "+str(rp_lim_val))
                    for ix in range(nx):
                     for iy in range(ny):
                      for iz in range(nz):
                          if(verbose):
                              print("Top find_clumps call at: "+str(ix)+"  "+str(iy)+"  "+str(iz))
                          find_clumps(ix,iy,iz,0)
                          
                    n_cl=len(cl_cells)
                    if(verbose):
                        print("Found "+str(n_cl)+" clumps")
                    N_clumps_entry[ie,xl,yl,tl,pl]=n_cl
                    fraction_volume=0.
                    singles=0
                    total_taken_cells=0
                    for q in range(len(cl_cells)):
                        volume=(cl_sizes[q][0][1]-cl_sizes[q][0][0]+1)*(cl_sizes[q][1][1]-cl_sizes[q][1][0]+1)*(cl_sizes[q][2][1]-cl_sizes[q][2][0]+1)
                        if(verbose):
                            print("volume: "+str(volume))
                        fraction_volume_new=cl_cells[q]/volume
                        total_taken_cells=total_taken_cells+cl_cells[q]
                        if(verbose):
                            print("number of cells: "+str(cl_cells[q]))
                        if(cl_cells[q]>1):
                          fraction_volume=fraction_volume+fraction_volume_new
                          if(verbose):
                              print("fraction volume: "+str(fraction_volume_new))
                        else:
                          singles=singles+1
                        if(verbose):
                            print("Volume occupancy of clump "+str(q)+": "+str(cl_cells[q]/volume))
                    N_clumps_single_cell_entry[ie,xl,yl,tl,pl]=singles
                    if(n_cl>0):
                        non_singles=n_cl-singles
                        if(verbose):
                            print("fraction volume total: "+str(fraction_volume))
                            print("Number of non singles: "+str(non_singles))
                        if(non_singles>0):
                            if(verbose):
                                print("Here")
                                print(type(fraction_volume))
                                print(type(non_singles))
                            N_clumps_percent_entry[ie,xl,yl,tl,pl]=fraction_volume
                            if(verbose):
                                print("Average volume occupancy: "+'{:7.5e}'.format(N_clumps_percent_entry[ie,xl,yl,tl,pl]))
                        N_clumps_single_cell_percent_entry[ie,xl,yl,tl,pl]=float(singles)/total_taken_cells
                        N_clumps_nonempty_entry[ie,xl,yl,tl,pl]=1
                    if(verbose):
                        print("Average volume occupancy: "+'{:7.5e}'.format(N_clumps_percent_entry[ie,xl,yl,tl,pl]))
                        print("total taken cells: "+str(total_taken_cells))
                        print("Single cells: "+str(singles))

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if(denominator[ix,iy,iz]>=p_min_for_non_emptyness):
                        Non_empty_cells=Non_empty_cells+1
                    else:
                        continue
                    for ie in range(len(elim)):
                        if(Tmunu[iT00,ix,iy,iz]<elim[ie]):
                            break
                        for xl in range(len(Xlim)):
                            if(X_ratio[ix,iy,iz]>Xlim[xl]):
                                break
                            for yl in range(len(Ylim)):
                                if(Y_ratio[ix,iy,iz]>Ylim[yl]):
                                    break
                                for tl in range(len(XT_lim)):
                                    if(XT_ratio[ix,iy,iz]>XT_lim[tl]):
                                         break
                                    for pl in range(len(RP_lim)):
                                         if(RP_entry[ix,iy,iz]>RP_lim[pl]):
                                              break
                                         else:
                                              N_cells_entry[ie,xl,yl,tl,pl]=N_cells_entry[ie,xl,yl,tl,pl]+1
        
        if(Non_empty_cells==0):
            print("Warning, it seems that all cells have zero energy density at time "+'{:7.5f}'.format(time))
            f_denominator=1
        else:
            f_denominator=Non_empty_cells
        if(n_i==0):
            T.append(Tmunu)
            en_per_part.append(energy_per_hadron)
            X.append(X_ratio)
            Y.append(Y_ratio)
            XT.append(XT_ratio)
            RP.append(RP_entry)
        else:
            T[index]=T[index]+Tmunu
            en_per_part[index]=en_per_part[index]+energy_per_hadron
            X[index]=X[index]+X_ratio
            Y[index]=Y[index]+Y_ratio
            XT[index]=XT[index]+XT_ratio
            RP[index]=RP[index]+RP_entry
        N_cells_event.append(N_cells_entry)
        N_cells_event_fraction.append(N_cells_entry/f_denominator)
        Non_empty_cells_event.append(Non_empty_cells)
        if(search_clumps):
            N_clumps_event.append(N_clumps_entry)
            N_clumps_single_cell_event.append(N_clumps_single_cell_entry)
            N_clumps_percent_event.append(N_clumps_percent_entry)
            N_clumps_single_cell_percent_event.append(N_clumps_single_cell_percent_entry)
            N_clumps_nonempty_event.append(N_clumps_nonempty_entry)

        index=index+1

        if(verbose):
            print("Done timestep: "+str(index)+", simulation time: "+str(time))

    #we check that we counted correctly the timesteps:
    if(len(tt)!=index):
        print("Error, I counted "+str(index)+" timesteps, but I have "+str(len(tt))+" entries in the list of timesteps...\nI quit.")
        sys.exit(2)
    if(n_i==0):
        nt=index
    N_cells_series.append(N_cells_event)
    N_cells_series_fraction.append(N_cells_event_fraction)
    Non_empty_cells_series.append(Non_empty_cells_event)
    if(search_clumps):
        N_clumps.append(N_clumps_event)
        N_clumps_nonempty.append(N_clumps_nonempty_event)
        N_clumps_single_cell.append(N_clumps_single_cell_event)
        N_clumps_percent.append(N_clumps_percent_event)
        N_clumps_single_cell_percent.append(N_clumps_single_cell_percent_event)

    fp_tmn.close()
    fp_jqbs.close()
    fp_vl.close()
    if(verbose):
        end_time = timer()
        print(i_tmn+", "+i_jqbs+", "+i_vl+" read in "+tf.format(end_time-start_time)+" seconds")
N_events=n_i+1

#we check that we counted correctly the number of events:
if(len(N_cells_series)!=N_events):
    print("Error, I counted "+str(N_events)+" events, but I have "+str(len(N_cells_series))+" sets of hydro candidate cells counting...\nI quit.")
    sys.exit(2)


N_cells_avg=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
N_cells_avg_fraction=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
Non_empty_cells_avg=np.zeros(nt,dtype=np.float64)
if(search_clumps):
    N_clumps_avg=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
    N_clumps_nonempty_tot=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
    N_clumps_single_cell_tot=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
    N_clumps_percent_tot=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
    N_clumps_single_cell_percent_tot=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)

if(verbose):
    print("Computing the averages.")
    start_time = timer()
for i in range(nt):
    T[i]=T[i]/N_events
    en_per_part[i]=en_per_part[i]/N_events
    X[i]=X[i]/N_events
    Y[i]=Y[i]/N_events
    XT[i]=XT[i]/N_events
    RP[i]=RP[i]/N_events
    for k in range(N_events):
        N_cells_avg[i,:,:,:,:,:]=N_cells_avg[i,:,:,:,:,:]+N_cells_series[k][i][:,:,:,:,:]
        N_cells_avg_fraction[i,:,:,:,:,:]=N_cells_avg_fraction[i,:,:,:,:,:]+N_cells_series_fraction[k][i][:,:,:,:,:]
        Non_empty_cells_avg[i]=Non_empty_cells_avg[i]+Non_empty_cells_series[k][i]

    if(search_clumps):    
        for k in range(N_events):
            N_clumps_avg[i,:,:,:,:,:]=N_clumps_avg[i,:,:,:,:,:]+N_clumps[k][i][:,:,:,:,:]
            N_clumps_nonempty_tot[i,:,:,:,:,:]=N_clumps_nonempty_tot[i,:,:,:,:,:]+N_clumps_nonempty[k][i][:,:,:,:,:]
            N_clumps_single_cell_tot[i,:,:,:,:,:]=N_clumps_single_cell_tot[i,:,:,:,:,:]+N_clumps_single_cell[k][i][:,:,:,:,:]
            N_clumps_percent_tot[i,:,:,:,:,:]=N_clumps_percent_tot[i,:,:,:,:,:]+N_clumps_percent[k][i][:,:,:,:,:]
            N_clumps_single_cell_percent_tot[i,:,:,:,:,:]=N_clumps_single_cell_percent_tot[i,:,:,:,:,:]+N_clumps_single_cell_percent[k][i][:,:,:,:,:]
    N_cells_avg[i,:,:,:,:,:]=N_cells_avg[i,:,:,:,:,:]/N_events
    N_cells_avg_fraction[i,:,:,:,:,:]=N_cells_avg_fraction[i,:,:,:,:,:]/N_events
    Non_empty_cells_avg[i]=Non_empty_cells_avg[i]/N_events
    if(search_clumps):
        N_clumps_avg[i,:,:,:,:,:]=N_clumps_avg[i,:,:,:,:,:]/N_events
    else: #we print these variables anyway, albeit in a way that minimize memory space usage
        N_clumps_avg=0
        N_clumps_nonempty_tot=0
        N_clumps_single_cell_tot=0
        N_clumps_percent_tot=0
        N_clumps_single_cell_percent_tot=0

if(verbose):
    end_time = timer()
    print("Done in "+tf.format(end_time-start_time)+" seconds")
#now we compute the standard deviations
N_cells_sd=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
N_cells_sd_fraction=np.zeros((nt,len(elim),len(Xlim),len(Ylim),len(XT_lim),len(RP_lim)),dtype=np.float64)
Non_empty_cells_sd=np.zeros(nt,dtype=np.float64)
if(verbose):
    print("Computing the standard deviations.")
    start_time = timer()
for i in range(nt):
    for k in range(N_events):
        N_cells_sd[i,:,:,:,:,:]=N_cells_sd[i,:,:,:,:,:]+(N_cells_avg[i,:,:,:,:,:]-N_cells_series[k][i][:,:,:,:,:])**2
        N_cells_sd_fraction[i,:,:,:,:,:]=N_cells_sd_fraction[i,:,:,:,:,:]+(N_cells_avg_fraction[i,:,:,:,:,:]-N_cells_series_fraction[k][i][:,:,:,:,:])**2
        Non_empty_cells_sd[i]=Non_empty_cells_sd[i]+(Non_empty_cells_avg[i]-Non_empty_cells_series[k][i])**2
    N_cells_sd[i,:,:,:,:,:]=np.sqrt(N_cells_sd[i,:,:,:,:,:]/(N_events-1))
    N_cells_sd_fraction[i,:,:,:,:,:]=np.sqrt(N_cells_sd_fraction[i,:,:,:,:,:]/(N_events-1))
    Non_empty_cells_sd[i]=np.sqrt(Non_empty_cells_sd[i]/(N_events-1))
if(verbose):
    end_time = timer()
    print("Done in "+tf.format(end_time-start_time)+" seconds")

# we transform the lists of 3D arrays into 4D arrays
# we do it one by one to save memory
T_arr=np.zeros((nt,10,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    T_arr[i,:,:,:]=T[i][:,:,:,:]
T=None
en_per_part_arr=np.zeros((nt,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    en_per_part_arr[i,:,:,:]=en_per_part[i][:,:,:]
en_per_part=None
X_arr=np.zeros((nt,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    X_arr[i,:,:,:]=X[i][:,:,:]
X=None
Y_arr=np.zeros((nt,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    Y_arr[i,:,:,:]=Y[i][:,:,:]
Y=None
XT_arr=np.zeros((nt,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    XT_arr[i,:,:,:]=XT[i][:,:,:]
XT=None
RP_arr=np.zeros((nt,nx,ny,nz),dtype=np.float64)
for i in range(nt):
    RP_arr[i,:,:,:]=RP[i][:,:,:]
RP=None

if(verbose):
    print("Writing the final results in "+outputfile)
    start_time = timer()

if print_only_energy_density_evolution:
    i0=3
    j0=4
    k0=1
    with open(outputfile,"w") as po:
      for i in range(nt):
        #po.write('{:4.2f}'.format(tt[i])+"    "+'{:14.9f}'.format(T[i][0,i0,j0,k0])+"    "+'{:14.9f}'.format(eps[i][i0,j0,k0])+"\n")
        po.write('{:4.2f}'.format(tt[i])+"    "+'{:14.9f}'.format(T[i][0,i0,j0,k0])+"\n")
else:
    with open(outputfile,"wb") as po:
      pickle.dump((lattice,tt,elim,Xlim,Ylim,XT_lim,RP_lim,N_events,en_per_part_arr,T_arr,X_arr,Y_arr,XT_arr,RP_arr,\
      N_cells_avg,N_cells_sd,Non_empty_cells_avg,Non_empty_cells_sd,\
      N_cells_avg_fraction,N_cells_sd_fraction,N_clumps_avg,N_clumps_percent_tot,N_clumps_single_cell_percent_tot,N_clumps_single_cell_tot,N_clumps_nonempty_tot),po)

if(verbose):
    end_time = timer()
    print("Done in "+tf.format(end_time-start_time)+" seconds")
    print("All done in "+tf.format(end_time-init_start)+" seconds")
