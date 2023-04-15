# version 0.1.1 - 17/12/2021 this program computes the averages and the standard deviations
#                          of energy densities, X_ebe and Y_ebe from thermodynamics.dat files
#                          In this basic version we use fixed parameters

import fileinput
import math
import numpy as np
import sys
import os
import matplotlib
import gzip
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

# warning, here the time arrays vary between different systems,
# there are offset, don't assume that the same time index between
# different arrays corresponds to the same time

#Tmunu indexes
T00=0
T01=1
T02=2
T03=3
T11=4
T12=5
T13=6
T22=7
T23=8
T33=9

# minimum pressure threshold [GeV/fm^3] to accept a point as part of the system
pmin=0.0001

# we leave this section commented for possible future use
# we get the name output directory
#N_input_files = len(sys.argv) - 1

#if N_input_files != 2:
#    print('Syntax: ./extract_data.py <input file> <output file>')
#    sys.exit(1)

#infile=sys.argv[1]
#if (not os.path.exists(infile)):
#    print("Input file "+infile+" does not exist.")
#    sys.exit(1)

#outfile = sys.argv[2]

outfile="time_evolution_comparison_w_stdev_at_center.pickle"
outfile_ascii="time_evolution_comparison_w_stdev_at_center.dat"

indirs=["res_Ag_Elab_1_58","res_Au_Elab_1_23_centr","res_Au_Ecm77"]
labels_ascii_file=["Ag+Ag, Elab=1.58AGeV, b=0-2.44fm ","Au+Au, Elab=1.23AGeV, b=0-3.3fm ","Au+Au, sqrt(s_NN)=7.7GeV, b=0-3.3fm "]
times_arr=[]
nt_arr=[]
nt_max=0
# we count the timesteps
for idir in indirs: 
  with open(idir+"/thermodynamics.dat","r") as ipf:
    times=[]
    count=False
    while(True):
        stuff=ipf.readline().split()
        if ((stuff[1]=="event") and (count==False)):
            count=True
            continue
        if ((stuff[1]=="event") and (count==True)):
            break
        if len(stuff)!=11:
            continue
        times.append(np.float64(stuff[0]))
    times_np=np.array(times)
    nt=times_np.size
    times_arr.append(times_np)
    nt_arr.append(nt)
    if (nt>nt_max):
         nt_max=nt

Tmunu=np.zeros((3,2,nt_max,10),dtype=np.float64)
accepted=np.zeros((3,nt_max),dtype=np.int32)
Xebe=np.zeros((3,2,nt_max),dtype=np.float64)
Yebe=np.zeros((3,2,nt_max),dtype=np.float64)
Xavg=np.zeros((3,nt_max),dtype=np.float64)
Yavg=np.zeros((3,nt_max),dtype=np.float64)

# first, we compute the averages
events_arr=[]
for ii,idir in enumerate(indirs): 
  infile=idir+"/thermodynamics.dat"
  with open(infile,"r") as ipf:
    #DBG print("HHH"+infile)
    n_events=0
    while(True):
        stuff=ipf.readline().split()
        if len(stuff) == 0:
            break
        if stuff[1]=="event":
            n_events=n_events+1
            #DBG print(str("EVENT: "+str(n_events-1)))
            for i in range(nt_arr[ii]):
                #DBG buff=ipf.readline().split()
                #DBG print(str(times_arr[ii][i])+"  "+str(buff[:]))
                Tmp=np.float64(ipf.readline().split())[1:]
                #DBG Tmp=np.float64(buff[1:])
                den=Tmp[T11]+Tmp[T22]+Tmp[T33]
                if den>=pmin:
                    X=(abs(Tmp[T11]-Tmp[T22])+abs(Tmp[T11]-Tmp[T33])+abs(Tmp[T33]-Tmp[T22]))/den
                    Y=3*(abs(Tmp[T13])+abs(Tmp[T23])+abs(Tmp[T12]))/den
                    accepted[ii,i]=accepted[ii,i]+1
                    #DBG if(i<3):
                        #DBG  print("AAAAAAAAAAAAAAAA")
                         #DBG print(str(n_events-1)+"  "+str(i)+"  "+str(den)+"  "+str(X)+"   "+str(Y))
                         #DBG print(str(Tmp[:]))
                         #DBG sys.exit(0)
                else:
                    X=0.
                    Y=0.
                    Tmp[:]=0.
                Xebe[ii,0,i]=Xebe[ii,0,i]+X
                Yebe[ii,0,i]=Yebe[ii,0,i]+Y
                Tmunu[ii,0,i,:]=Tmunu[ii,0,i,:]+Tmp[:]

  print(infile+" read, events: "+str(n_events))
  events_arr.append(n_events)
  Tmunu[ii,0,:,:]=Tmunu[ii,0,:,:]/n_events
  for i in range(nt_arr[ii]):
      if accepted[ii,i]>0 :
          Xebe[ii,0,i]=Xebe[ii,0,i]/accepted[ii,i]
          Yebe[ii,0,i]=Yebe[ii,0,i]/accepted[ii,i]

# now we compute the standard deviations
for ii,idir in enumerate(indirs): 
  infile=idir+"/thermodynamics.dat"
  with open(infile,"r") as ipf:
    while(True):
        stuff=ipf.readline().split()
        if len(stuff) == 0:
            break
        if stuff[1]=="event":
            for i in range(nt_arr[ii]):
                Tmp=np.float64(ipf.readline().split())[1:]
                den=Tmp[T11]+Tmp[T22]+Tmp[T33]
                if den>=pmin:
                    X2=(Xebe[ii,0,i]-(abs(Tmp[T11]-Tmp[T22])+abs(Tmp[T11]-Tmp[T33])+abs(Tmp[T33]-Tmp[T22]))/den)**2
                    Y2=(Yebe[ii,0,i]-3*(abs(Tmp[T13])+abs(Tmp[T23])+abs(Tmp[T12]))/den)**2
                else:
                    X2=0.
                    Y2=0.
                Xebe[ii,1,i]=Xebe[ii,1,i]+X2
                Yebe[ii,1,i]=Yebe[ii,1,i]+Y2
                Tmunu[ii,1,i,:]=Tmunu[ii,1,i,:]+(Tmp[:]-Tmunu[ii,0,i,:])**2

  for i in range(nt_arr[ii]):
      Tmunu[ii,1,i,:]=np.sqrt(Tmunu[ii,1,i,:]/events_arr[ii])
      if accepted[ii,i]>0 :
          Xebe[ii,1,i]=np.sqrt(Xebe[ii,1,i]/accepted[ii,i])
          Yebe[ii,1,i]=np.sqrt(Yebe[ii,1,i]/accepted[ii,i])

      den=Tmunu[ii,0,i,T11]+Tmunu[ii,0,i,T22]+Tmunu[ii,0,i,T33]
      if den>=pmin:
         Xavg[ii,i]=(abs(Tmunu[ii,0,i,T11]-Tmunu[ii,0,i,T22])+abs(Tmunu[ii,0,i,T11]-Tmunu[ii,0,i,T33])+abs(Tmunu[ii,0,i,T33]-Tmunu[ii,0,i,T22]))/den
         Yavg[ii,i]=3*(abs(Tmunu[ii,0,i,T13])+abs(Tmunu[ii,0,i,T23])+abs(Tmunu[ii,0,i,T12]))/den

with open(outfile,"wb") as outf:
    pickle.dump((times_arr,events_arr,Tmunu,Xebe,Yebe,Xavg,Yavg),outf)

sp="    "
with open(outfile_ascii,"w") as outf:
    outf.write("# for each system, after printing its basic properties (ion types, collision energy, centrality)\n")
    outf.write("# we print, for its central point: \n")
    outf.write("# column 1: time [fm]\n")
    outf.write("# column 2-11: Tmunu components 2=T00, 3=T01, 4=T02, 5=T03, 6=T11, 7=T12, 8=T13, 9=T22, 10=T23, 11=T33\n")
    outf.write("# column 12-21: standard deviations of Tmunu components 12=T00, 13=T01, 14=T02, 15=T03, 16=T11, 17=T12, 18=T13, 19=T22, 20=T23, 21=T33\n")
    outf.write("# colum 22: X_ebe, column 23: standard deviation of X_ebe\n")
    outf.write("# colum 24: Y_ebe, column 25: standard deviation of Y_ebe\n")
    outf.write("# column 26: X, column 27: Y\n")
    for ii,ilab in enumerate(labels_ascii_file):
         outf.write("\n# system: "+ilab+"\n")
         for i in range(nt_arr[ii]):
             outf.write('{:3.1f}'.format(times_arr[ii][i])+sp)
             for k in range(10):
                 outf.write('{:9.5e}'.format(Tmunu[ii,0,i,k])+sp)
             for k in range(10):
                 outf.write('{:9.5e}'.format(Tmunu[ii,1,i,k])+sp)
             outf.write('{:9.5e}'.format(Xebe[ii,0,i])+sp+'{:9.5e}'.format(Xebe[ii,1,i])+sp)
             outf.write('{:9.5e}'.format(Yebe[ii,0,i])+sp+'{:9.5e}'.format(Yebe[ii,1,i])+sp)
             outf.write('{:9.5e}'.format(Xavg[ii,i])+sp+'{:9.5e}'.format(Yavg[ii,i])+"\n")
