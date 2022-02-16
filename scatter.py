#  ---------------------------------------------------------------------------------------------  #
#
#   The program plots the sensitivities for all the values of K found in Results/
#
#  ---------------------------------------------------------------------------------------------  #

import sys
import numpy as np
from matplotlib import pyplot as plt
from glob import glob

T = float( sys.argv[1] )
deltat = float( sys.argv[2] )

L = int(T/deltat)

def filename_func():
    return "Results/T%.4f_dt%.4f_K?.????.txt" % (T, deltat)

files = glob(filename_func())
Ks = np.array([float(f[-10:-4]) for f in files])

ordr = np.argsort(Ks)

res = np.zeros((len(Ks),5)) # K <pulses> std(pulses) <1/eta> std(1/eta) 

for k in range(len(files)):
    
    data = np.loadtxt(files[ordr[k]]).T
    plt.plot(data[0], data[1], 'o', ms=2, label="K=%.1e"%Ks[ordr[k]])
    
    res[k,0] = Ks[ordr[k]]
    res[k,1] = np.average(data[0])
    res[k,2] = np.std(data[0])
    res[k,3] = np.average(data[1])
    res[k,4] = np.std(data[1])

np.savetxt("varK_T%.4f.txt" % T, res, header="K <pulses> std(pulses) <1/eta> std(1/eta)")

plt.xlabel("# of pulses")
plt.ylabel(r"1/$\eta$")

#plt.title(r"$T=%.2f \mu$s, monochromatic signal" % T)
plt.title(' '.join(sys.argv[1:]))

plt.legend()
plt.show()