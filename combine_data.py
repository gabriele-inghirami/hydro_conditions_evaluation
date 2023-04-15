# combine_data.py - version 0.8.4 - 14/10/2021

# it works with analize_data v 0.8.x

import fileinput
import math
import numpy as np
import sys
import os
import pickle
from timeit import default_timer as timer

#if False it prints only error messages, if True it writes what it is doing at the moment and the intermediate results 
verbose=True

#function that combines two variances
def combvar(m1, m2, v1, v2, n1, n2):
    #m1 and m2 are the two means, v1 and v2 the two variances, n1 and n2 the number of elements in the two sets
    return (v1*(n1-1)+v2*(n2-1))/(n1+n2-1) + n1*n2*((m1-m2)**2)/((n1+n2)*(n1+n2-1))


if(verbose):
    init_start=timer()

#we parse the command line arguments
N_input_args=len(sys.argv)-1

if(N_input_args<2):
   print ('Syntax: ./combine_data.py <file data 1> [data 2] ... <outputfile>')
   print ("file data 1,2,3...N are the pickled files produced by analize_data.py")
   print ("outputfile is obviously the name of the output file with the results of the postprocessing")
   sys.exit(1)

#we get the name of input and output files
inputfiles=sys.argv[1:N_input_args]
n_input=len(inputfiles)
outputfile=sys.argv[N_input_args]

for n_i, infile in enumerate(inputfiles):
    if(infile[-3:]==".gz"):
       if(verbose):
           print("Opening gzipped file "+infile)
       pi=gzip.open(ifile,"rb")
    else:
        if(verbose):
             print("Opening file "+infile)
        pi=open(infile,"rb")

    indata=pickle.load(pi)
    pi.close()
    lattice,tt,elim,Xlim,Ylim,XT_lim,RP_lim,N_events,en_per_part,Tmunu,X,Y,XT,RP,N_cells_avg,N_cells_sd,\
    Non_empty_cells_avg,Non_empty_cells_sd,N_cells_avg_fraction,N_cells_sd_fraction,N_clumps_avg,N_clumps_percent_tot,N_clumps_single_cell_percent_tot,\
    N_clumps_single_cell_tot,N_clumps_nonempty_tot = indata[:]


    if(isinstance(N_clumps_avg,np.ndarray)):
        clumps_on=True
    else:
        clumps_on=False

    if(n_i==0):
        lattice_ref=lattice
        tt_ref=tt
        nt=len(tt)
        elim_ref=elim
        Xlim_ref=Xlim
        Ylim_ref=Ylim
        XT_lim_ref=XT_lim
        RP_lim_ref=RP_lim
        N_events_total=N_events
        Tmunu_total=Tmunu*N_events
        en_per_part_total=en_per_part*N_events
        X_total=X*N_events
        Y_total=Y*N_events
        XT_total=XT*N_events
        RP_total=RP*N_events
        N_cells_total=N_cells_avg*N_events
        variance_N_cells_total=N_cells_sd**2
        Non_empty_cells_total=Non_empty_cells_avg*N_events
        variance_total_non_empty_cells=Non_empty_cells_sd**2
        N_cells_fraction_total=N_cells_avg_fraction*N_events
        variance_total_fraction=N_cells_sd_fraction**2
        if(clumps_on):
            N_clumps_total=N_clumps_avg*N_events
            N_clumps_single_cell_total=N_clumps_single_cell_tot
            N_clumps_nonempty_total=N_clumps_nonempty_tot
            N_clumps_percent_total=N_clumps_percent_tot
            N_clumps_single_cell_percent_total=N_clumps_single_cell_percent_tot
    else:
        if(not np.array_equal(lattice_ref["dimensions"],lattice["dimensions"])):
            print("Error in file "+infile+": different lattice dimensions. Until now: "+str(lattice_ref["dimensions"])+", this time: "+str(lattice["dimensions"])+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice_ref["spacing"],lattice["spacing"])):
            print("Error in file "+infile+": different lattice spacing. Until now: "+str(lattice_ref["spacing"])+", this time: "+str(lattice["spacing"])+".\nI quit.")
            sys.exit(2)
        if(not np.array_equal(lattice_ref["origin"],lattice["origin"])):
            print("Error in file "+infile+": different lattice origin. Until now: "+str(lattice_ref["origin"])+", this time: "+str(lattice["origin"])+".\nI quit.")
            sys.exit(2)
        if(elim!=elim_ref):
            print("Error in file "+infile+": different energy density tresholds. Until now: "+str(elim_ref)+", this time: "+str(elim)+".\nI quit.")
            sys.exit(2)
        if(Xlim!=Xlim_ref):
            print("Error in file "+infile+": different energy density tresholds. Until now: "+str(Xlim_ref)+", this time: "+str(Xlim)+".\nI quit.")
            sys.exit(2)
        if(Ylim!=Ylim_ref):
            print("Error in file "+infile+": different energy density tresholds. Until now: "+str(Ylim_ref)+", this time: "+str(Ylim)+".\nI quit.")
            sys.exit(2)
        if(XT_lim!=XT_lim_ref):
            print("Error in file "+infile+": different energy density tresholds. Until now: "+str(XT_lim_ref)+", this time: "+str(XT_lim)+".\nI quit.")
            sys.exit(2)
        if(RP_lim!=RP_lim_ref):
            print("Error in file "+infile+": different energy density tresholds. Until now: "+str(RP_lim_ref)+", this time: "+str(RP_lim)+".\nI quit.")
            sys.exit(2)
        N_events_old=N_events_total
        N_events_total=N_events_old+N_events
        Tmunu_total=Tmunu_total+Tmunu*N_events
        en_per_part_total=en_per_part_total+en_per_part*N_events
        X_total=X_total+X*N_events
        Y_total=Y_total+Y*N_events
        XT_total=XT_total+XT*N_events
        RP_total=RP_total+RP*N_events
        N_cells_avg_old=N_cells_total/N_events_old
        N_cells_total=N_cells_total+N_cells_avg*N_events
        N_cells_avg_fraction_old=N_cells_fraction_total/N_events_old
        N_cells_fraction_total=N_cells_fraction_total+N_cells_avg_fraction*N_events
        Non_empty_cells_avg_old=Non_empty_cells_total/N_events_old
        Non_empty_cells_total=Non_empty_cells_total+Non_empty_cells_avg*N_events

        variance_N_cells_total_old=variance_N_cells_total
        variance_N_cells_total_new=N_cells_sd**2
        variance_N_cells_total=combvar(N_cells_avg_old, N_cells_avg, variance_N_cells_total_old, variance_N_cells_total_new, N_events_old, N_events)

        variance_total_fraction_old=variance_total_fraction
        variance_total_fraction_new=N_cells_sd_fraction**2
        variance_total_fraction=combvar(N_cells_avg_fraction_old, N_cells_avg_fraction, variance_total_fraction_old, variance_total_fraction_new, N_events_old, N_events)

        variance_total_non_empty_cells_old=variance_total_non_empty_cells
        variance_total_non_empty_cells_new=Non_empty_cells_sd**2
        variance_total_non_empty_cells=combvar(Non_empty_cells_avg_old, Non_empty_cells_avg, variance_total_non_empty_cells_old, variance_total_non_empty_cells_new, N_events_old, N_events)

        if(clumps_on):
            N_clumps_total=N_clumps_total+N_clumps_avg*N_events
            N_clumps_single_cell_total=N_clumps_single_cell_total+N_clumps_single_cell_tot
            N_clumps_nonempty_total=N_clumps_nonempty_total+N_clumps_nonempty_tot
            N_clumps_percent_total=N_clumps_percent_total+N_clumps_percent_tot
            N_clumps_single_cell_percent_total=N_clumps_single_cell_percent_total+N_clumps_single_cell_percent_tot

N_events=N_events_total
Tmunu=Tmunu_total/N_events
en_per_part=en_per_part_total/N_events
X=X_total/N_events
Y=Y_total/N_events
XT=XT_total/N_events
RP=RP_total/N_events
N_cells_avg=N_cells_total/N_events
N_cells_sd=np.sqrt(variance_N_cells_total)
Non_empty_cells_avg=Non_empty_cells_total/N_events
Non_empty_cells_sd=np.sqrt(variance_total_non_empty_cells)
N_cells_avg_fraction=N_cells_fraction_total/N_events
N_cells_sd_fraction=np.sqrt(variance_total_fraction)

if(clumps_on):
    N_clumps_avg=N_clumps_total/N_events
else: #we print these variables anyway, albeit in a way that minimize memory space usage
    N_clumps_avg=0
    N_clumps_nonempty_total=0
    N_clumps_single_cell_total=0
    N_clumps_percent_total=0
    N_clumps_single_cell_percent_total=0


with open(outputfile,"wb") as po:
    pickle.dump((lattice,tt,elim,Xlim,Ylim,XT_lim,RP_lim,N_events,en_per_part,Tmunu,X,Y,XT,RP,N_cells_avg,N_cells_sd,\
    Non_empty_cells_avg,Non_empty_cells_sd,N_cells_avg_fraction,N_cells_sd_fraction,N_clumps_avg,N_clumps_percent_total,N_clumps_single_cell_percent_total,N_clumps_single_cell_total,N_clumps_nonempty_total),po)
