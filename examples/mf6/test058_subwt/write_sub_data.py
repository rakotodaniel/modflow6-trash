import numpy as np
import flopy

ml = flopy.modflow.Modflow.load("subwt_mf2005.nam",model_ws="mf2005")

# calculate frac
lay_thick = ml.dis.thickness.array
compress_thick = [45.,57.,50.,90.]
for i,ct in enumerate(compress_thick):
    lay_thick[i,:,:] = ct/lay_thick[i,:,:]

# get node numbers
kijs,nns = [],[]
nrc = ml.nrow * ml.ncol
for k in range(ml.nlay):
    for i in range(ml.nrow):
        kijs.extend([(k,i,j) for j in np.arange(ml.ncol)])
        nns.extend([int(((k) * nrc) + ((i) * ml.ncol) + j) for j in np.arange(ml.ncol)])

cc = 0.25
cr = 0.01
void = 0.82
kv = -999
sgm = 1.7
sgs = 2.0
ini_stress = 15.0
delay_flag = 0
idomain = np.loadtxt("idomain.dat")

with open("subwt_mf6.sub.data",'w') as f:
    f.write("#l r c initial_stress frac cc_sse cr_ssv void isdelay boundname\n")
    for kij,nn in zip(kijs,nns):
        if idomain[kij[1],kij[2]] == 0:
            continue
        f.write("{0} {1} {2}".format(kij[0]+1,kij[1]+1,kij[2]+1))
        f.write(" {0} {1:.3f} {2} {3} {4} {5} {6} nsystm{7} \n".format(ini_stress,lay_thick[kij],
            cc,cr,void,kv,delay_flag,kij[0]))
        

