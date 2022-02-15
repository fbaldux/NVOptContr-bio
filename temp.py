#  ---------------------------------------------------------------------------------------------  #
#
#   Given the final time and tone, the program plots the sensitivities for all the
#   values of K found in Results/
#
#  ---------------------------------------------------------------------------------------------  #

import sys
import numpy as np
from matplotlib import pyplot as plt

T = 100
deltat = 0.02


data = np.loadtxt("varK_T%.4f.txt" % T).T

which = data[0]>0
plt.errorbar(data[1,which], data[3,which], xerr=data[2,which], yerr=data[4,which], fmt='o', ms=3, c='black', label="vanilla")

which = data[0]==0
plt.errorbar(data[1,which], data[3,which], xerr=data[2,which], yerr=data[4,which], fmt='s', ms=3, c='red', label="spherical")


plt.xlabel("# of pulses")
plt.ylabel(r"1/$\eta$")

plt.title(r"$T=%.2f \mu$s, trying to reduce $\pi$ pulses" % T)

plt.legend()
plt.savefig("plot.pdf", bbox_inches='tight')
plt.show()