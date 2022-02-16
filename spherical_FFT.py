#  ---------------------------------------------------------------------------------------------  #
#   
#   The program finds the configuration of continuous spins s[i] that minimezes the cost function
#       
#       H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - sum_i s[i]**2
#
#   - The variables `J[i,j]` and `h[i]` are loaded from Init/
#   - From the continuous spins are generated Ising spins s_Ising[i] = sign(s[i]), that are 
#     then saved to Configurations/, with the # of pulses and 1/eta in the header.
#
#  ---------------------------------------------------------------------------------------------  #


import sys
import numpy as np
from scipy.linalg import toeplitz
from scipy.fft import fft,ifft,fftshift
from scipy.optimize import brentq
from matplotlib import pyplot as plt
from time import time


Tfin = float( sys.argv[1] )
Delta_t = float( sys.argv[2] )
alpha = float( sys.argv[3] )

N = int( Tfin / Delta_t )


#  ------------------------------------------  import  -----------------------------------------  #

J = np.loadtxt("Init/J_T%.4f_dt%.4f_a%.4f.txt" % (Tfin,Delta_t,alpha))
Jmat = toeplitz(J)

# to have decay both for positive and negative indices
#J[N//2+1:] = np.flip(J[1:N//2])

h = np.loadtxt("Init/h_T%.4f_dt%.4f.txt" % (Tfin,Delta_t))


#  ------------------------------------  Fourier transform  ------------------------------------  #

start = time()

hF = fft(h, norm='ortho')
JF = fft(J, norm='ortho')

#plt.plot(np.arange(len(h)), h)
"""plt.plot(np.arange(len(hF)), np.abs(fftshift(hF)))
plt.yscale("log")
plt.show()
exit(0)"""


#  ---------------------------------  functions for S.P. eqs.  ---------------------------------  #

# l.h.s. of the equation
def lhsTemp(lamda):
    return np.sum( np.abs( hF / (JF + lamda) )**2 )
lhs = np.vectorize(lhsTemp)

# r.h.s. of the equation
def rhsTemp(lamda):
    return N * np.abs( np.sum( np.abs(hF)**2 / (JF + lamda) ) )
rhs = np.vectorize(rhsTemp)

# combined equation
def equation(lamda):
    return lhs(lamda) - rhs(lamda)
    

#  -------------------------------------  solve S.P. eqs.  -------------------------------------  #

# find interval [a,b] in which equation changes sign
a = 0
b = 0.1
sign0 = np.sign(lhs(a)-rhs(a))

while np.sign(lhs(b)-rhs(b)) == sign0:
    b += 0.1

# solve equation in [a,b]
lamda = brentq(equation, a, b)


#  -------------------------------  transform back to real space  ------------------------------  #

# non-ortho
#C = np.sqrt( N * np.sum( np.abs(hF)**2/(JF+lamda) ) )
#sF = N*hF / (C*(JF+lamda))

# ortho
C = np.sqrt( np.sum( np.abs(hF)**2/(JF+lamda) ) )
sF = hF / ( C*(JF+lamda) )

s = ifft(sF, norm='ortho').real

s_Ising = np.sign(s).astype(np.int_)

"""
ft = open("times.txt", 'a')
ft.write("%f " % (time()-start))
ft.close()
"""

#  -------------------------------------------  plot  ------------------------------------------  #
"""
#s_ann = np.loadtxt("Configurations/s_T%.4f_dt%.4f_t%d_h%d_K0.0050_r0.txt" % (Tfin,Delta_t,tone,harmonic))

plt.plot(np.arange(N), s, '-', c='black', label=r"$s_i \in \mathbb{R}$")
plt.plot(np.arange(N), s_Ising, '--', c='firebrick', label=r"$s_i = \pm 1$")
#plt.plot(np.arange(N), s_ann, ':', c='darkgreen', label=r"sim. anneal.")
plt.plot(np.arange(N), h*np.sqrt(N/np.sum(h**2)), ':', c='darkgreen', label=r"$h(t)$")

plt.xlabel(r"$i$")
#plt.ylabel(r"$s_i$")

plt.title(r"solution for $N=%d$ spins, $T_f = %.1f$" % (N,Tfin))

plt.legend()
plt.show()
#exit(0)
"""
#  ------------------------------------  sensitivity & co.  ------------------------------------  #

def energy(s):
    return 0.5*np.einsum("a,ab,b", s, Jmat, s) - np.log( np.abs( np.dot(h,s) ) )

def domain_walls(s):
    #return (N - np.dot(s[:-1],s[1:]) - 1) // 2
    return np.sum( np.abs(np.diff(s)) ) // 2

def etaInv(epsilon):
    return 1. / np.exp( epsilon - np.log(28025) - 0.5*np.log(Tfin*1e-6) )


#print("phi", np.dot(h,s_Ising))
#print("chi", 0.5*np.einsum("a,ab,b", s_Ising, Jmat, s_Ising))


#  ---------------------------------------  save to file  --------------------------------------  #

filename = "Configurations/sSpher_T%.4f_dt%.4f_a%.4f.txt" % (Tfin,Delta_t,alpha)
head = "pulses=%d, 1/eta=%f" % (domain_walls(s_Ising), etaInv(energy(s_Ising)))
np.savetxt(filename, s_Ising, header=head, fmt='%d')




#print(Delta_t, etaInv(energy(s)), etaInv(energy(s_Ising)))
#print("python:", Delta_t, energy(s_Ising))
#print( "%f %e %e %e %e" % (Delta_t, energy(s), energy(s_Ising), etaInv(energy(s)), etaInv(energy(s_Ising))) )


















