# version 0.1.1 - 13/10/2021 this program averages the Tmunu computed in VTK format

import fileinput
import math
import numpy as np
import sys
import os

gs=29
g2=gs*gs
g3=g2*gs
origin='{:3.1f}'.format(-gs/2.) #"-1.5"
data=np.zeros((10,g2,gs),dtype=np.float64)
gs_string='{:3d}'.format(gs)
g3_string='{:5d}'.format(g3)

# we get the name output directory
N_input_files = len(sys.argv) - 2

if N_input_files < 1:
    print('Syntax: ./extract_data.py <output file> <input file(s)>')
    sys.exit(1)

outfile = sys.argv[1]

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

strings=[]
for i in range(N_input_files):
  with open(sys.argv[i+2],"r") as ipf:
        for k in range(0,8):
            strings.append(ipf.readline())
        for k in range(0,10):
            strings.append(ipf.readline())
            strings.append(ipf.readline())
            for j in range(g2):
                data[k,j,:]=data[k,j,:]+np.float64(ipf.readline().split())
        print("File "+sys.argv[i+2] +" read (#"+str(i)+")")


data=data/N_input_files
strings.reverse()

with open(outfile,"w") as ipf:
        for k in range(0,8):
            ipf.write(strings.pop())
        for k in range(0,10):
            ipf.write(strings.pop())
            ipf.write(strings.pop())
            for j in range(g2):
                for l in range(gs):
                    ipf.write('{:14.10f}'.format(data[k,j,l]))
                ipf.write("\n")


den=data[iT11,:,:]+data[iT22,:,:]+data[iT33,:,:]
num=abs(data[iT11,:,:]-data[iT22,:,:])+abs(data[iT11,:,:]-data[iT33,:,:])+abs(data[iT22,:,:]-data[iT33,:,:])
X_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den!=0)
num=3*(abs(data[iT12,:,:])+abs(data[iT13,:,:])+abs(data[iT23,:,:]))
Y_avg=np.divide(num,den,out=np.zeros(num.shape, dtype=np.float64), where=den!=0)

def write_output(filename,quantity):
  with open(filename,"w") as ipf:
     ipf.write("# vtk DataFile Version 2.0\n")
     ipf.write("hadron_tmn_landau\n")
     ipf.write("ASCII\n")
     ipf.write("DATASET STRUCTURED_POINTS\n")
     ipf.write("DIMENSIONS"+gs_string+gs_string+gs_string+"\n")
     ipf.write("SPACING 1 1 1\n")
     ipf.write("ORIGIN "+origin+" "+origin+" "+origin+"\n")
     ipf.write("POINT_DATA "+g3_string+"\n")
     ipf.write("SCALARS X_p_an double 1\n")
     ipf.write("LOOKUP_TABLE default\n")
     for j in range(g2):
         for l in range(gs):
             ipf.write('{:14.10f}'.format(quantity[j,l]))
         ipf.write("\n")

write_output("X.vtk",X_avg)
write_output("Y.vtk",Y_avg)
