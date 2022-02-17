import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

#Tfin = np.array((1,2,3,4,5,10,20,30,40,50,100))
Tfin = np.arange(2,17,2)
Delta_t = 0.01

#alpha = np.array((0.5,1,1.5,))
#alpha = np.array((2,))
a = 1

amplNoise = np.array((0.01,0.1,0.5,1))
#A = 1

cols = cm.get_cmap('magma', 5)

c = 0
#for ia in range(len(alpha)):
#    a = alpha[ia]
for iA in range(len(amplNoise)):
    A = amplNoise[iA]

    first = True
    for iT in range(len(Tfin)):
        T = Tfin[iT]
    
        filename = "Results/T%.4f_dt%.4f_a%.4f_A%.2e_K%.4f.txt" % (T,Delta_t,a,A,0)
        data = np.loadtxt(filename).T

        if first:
            #lab = r"$\alpha$=%.1f"%a
            lab = r"$A$=%.4f" % A
            plt.plot(T, np.average(data[1]), '.', c=cols(c), label=lab)
            first = False
        else:
            plt.plot(T, np.average(data[1]), '.', c=cols(c))
            
    
    c += 1

plt.xlabel(r"$T$ [$\mu$s]")
plt.ylabel(r"$1/\eta$ [Hz$^{1/2}$/$\mu$T]")
#plt.ylabel(r"# of pulses")

#plt.yscale("log")

#plt.ylim((0,40))

plt.title(r"annealing from spherical: $\Delta t$=%.2f" % Delta_t)

plt.legend()
plt.show()