import numpy as np
from matplotlib import pyplot as plt


#Tfin = np.array((1,2,3,4,5,10,20,30,40,50,100))
Tfin = np.arange(10,310,10)
Delta_t = 0.02


for iT in range(len(Tfin)):
    T = Tfin[iT]
    
    filename = "Results/T%.4f_dt%.4f_K%.4f.txt" % (T,Delta_t,0)
    data = np.loadtxt(filename).T

    plt.plot(T, data[0], '.', c='black')

plt.xlabel(r"$T$ [$\mu$s]")
#plt.ylabel(r"$1/\eta$")
plt.ylabel(r"# of pulses")

plt.title("annealing from spherical")

plt.show()