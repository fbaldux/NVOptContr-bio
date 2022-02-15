#  ---------------------------------------------------------------------------------------------  #
#
#   Given the final time and tone, the program plots the sensitivities for all the
#   values of K found in Results/
#
#  ---------------------------------------------------------------------------------------------  #

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

Tfin = 1000
Delta_t = 0.5

Ks = np.array([5.00e-03,1.48e-02,4.39e-02,1.30e-01,3.85e-01,1.14e+00,3.38e+00,1.00e+01])


#plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 16})
#plt.rcParams.update({"text.usetex": True})
#plt.rcParams["figure.figsize"] = [5,5]
fig, ax = plt.subplots()

cols = cm.get_cmap("inferno_r", len(Ks)+1)


def fexp(f):
    return int(np.floor(np.log10(abs(f)))) if f != 0 else 0

def fman(f):
    return f/10**fexp(f)
    

###################  ONLY ANNEALING  ###################

for iK in range(len(Ks)):
    K = Ks[iK]
    
    filename = "Results/T%.4f_dt%.4f_K%.4f.txt" % (Tfin,Delta_t,K)
    data = np.loadtxt(filename).T
    #ax.plot(data[0], data[1], 'o', ms=2, label="K=%.1e"%K, c=cols(iK))
    
    av_x = np.average(data[0])
    std_x = np.std(data[0])
    av_y = np.average(data[1])
    std_y = np.std(data[1])
    
    lab = r"$%.2f$$\times$$10^{%d}$" % (fman(K), fexp(K))
    ax.errorbar(av_x, av_y, xerr=std_x, yerr=std_y, marker='o', ms=4, label=lab, c=cols(iK+1))


###################  EXACT -> ANNEALING  ###################

ax.plot(0, 94.220404, marker='s', ms=4, c='darkgreen')

ax.text(10,94, "spherical", c='darkgreen')

###################  PLOT  ###################

    
ax.set_xlabel(r"$\pi$-pulse number")
ax.set_ylabel(r"$1/\eta$ [Hz$^{1/2}$/$\mu$T]")


#ax.set_xscale('log')

ax.legend(ncol=2, loc='lower right', title=r"$K$")
plt.show()