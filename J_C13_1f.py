#  ---------------------------------------------------------------------------------------------  #
#
#   Program that does the same of the .cpp version, just for checks.
#
#  ---------------------------------------------------------------------------------------------  #


import sys
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

Tfin = float( sys.argv[1] )
Delta_t = float( sys.argv[2] )
alpha = float( sys.argv[3] )
amplNoise = float( sys.argv[4] )
N = int( Tfin / Delta_t )


#  ----------------------------------------  functions  ----------------------------------------  #

# noise power spectrum: peaked contribution
def S1(omega):
    return 0.52 * np.exp(-0.5 * (omega-2.7156)*(omega-2.7156) * 1436)

# noise power spectrum: flat contribution
def S2(omega):
    return 0.00119;

# noise power spectrum: 1/f^alpha contribution
def S3(omega):
    return amplNoise * omega**(-alpha)


# full integrand: x = omega*Delta_t
def integrand(x,S,k):
    return S(x/Delta_t) * (1-np.cos(x)) * np.cos(k*x) / x**2


def J_int(S,k,xMin,xMax):
    I = quad(integrand, xMin, xMax, args=(S,k))[0]
    return 4*Delta_t/np.pi * I


#  -------------------------------------  computation of J  ------------------------------------  #

Js = np.zeros(N)

# numerical integration
for k in range(N):
    #Js[k] = J_int(S1, k, 2.4517*Delta_t, 2.9795*Delta_t)
    #Js[k] += J_int(S2, k, 0.001*Delta_t, 6.5*Delta_t)
    
    if amplNoise>0:
        Js[k] += J_int(S3, k, 0.001*Delta_t, 35*Delta_t)


#  --------------------------------------  saving to file  -------------------------------------  #
"""
filename = "Init/J_T%.4f_dt%.4f_a%.4f_A%.2e.txt" % (Tfin,Delta_t,alpha,amplNoise)
np.savetxt(filename, Js)
"""

#  ------------------------------------  comparison w/ C++  ------------------------------------  #


filename = "Init/J_T%.4f_dt%.4f_a%.4f_A%.2e.txt" % (Tfin,Delta_t,alpha,amplNoise)
J2 = np.loadtxt(filename)

plt.plot(np.arange(N), Js, 'o', ms=4, label="buono")
plt.plot(np.arange(N), J2, '+', ms=4, label="C++")

#plt.yscale("log")

plt.legend()
plt.show()
