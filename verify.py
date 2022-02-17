import sys
import numpy as np
from scipy.linalg import toeplitz
#from scipy.fft import fft,ifft,fftshift,fftfreq
from matplotlib import pyplot as plt

Tfin = float( sys.argv[1] )
Delta_t = float( sys.argv[2] )
alpha = float( sys.argv[3] )
amplNoise = float( sys.argv[4] )
K = float( sys.argv[5] )
N = int( Tfin / Delta_t )


#  ------------------------------------------  noise  ------------------------------------------  #

filename = "Init/J_T%.4f_dt%.4f_a%.4f_A%.2e.txt" % (Tfin,Delta_t,alpha,amplNoise)
J = np.loadtxt(filename)
Jmat = toeplitz(J)

plt.plot(np.arange(N), J, '.')

#  ------------------------------------------  signal  -----------------------------------------  #

filename = "Init/h_T%.4f_dt%.4f.txt" % (Tfin,Delta_t)
h = np.loadtxt(filename)

plt.plot(np.arange(N), h, '.')

plt.show()

#  ------------------------------------------  filter  -----------------------------------------  #

filename = "Configurations/s_T%.4f_dt%.4f_a%.4f_A%.2e_K%.4f_r%d.txt" % (Tfin,Delta_t,alpha,amplNoise,K,0)
s = np.loadtxt(filename)


#  ------------------------------------  sensitivity & co.  ------------------------------------  #

def energy(s):
    Jterm = 0.5*np.einsum("a,ab,b", s, Jmat, s)
    hterm = np.abs( np.dot(h,s) )
    print("Jterm =", Jterm)
    print("hterm =", hterm, "-> -log(hterm):", -np.log(hterm))
    return Jterm - np.log(hterm)

def domain_walls(s):
    return np.sum( np.abs(np.diff(s)) ) // 2

def etaInv(epsilon):
    return 1. / np.exp( epsilon - np.log(28025) - 0.5*np.log(Tfin*1e-6) )


print("N_pulses=%d; 1/eta=%f" % (domain_walls(s), etaInv(energy(s))))