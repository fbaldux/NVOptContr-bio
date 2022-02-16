import sys
import numpy as np
from scipy.fft import fft,ifft,fftshift,fftfreq
from matplotlib import pyplot as plt

Tfin = float( sys.argv[1] )
Delta_t = float( sys.argv[2] )
alpha = float( sys.argv[3] )

N = int( Tfin / Delta_t )

A = 1


#  ----------------------------------------  utilities  ----------------------------------------  #

def gaussian(x,sigma2):
    return np.exp(-0.5*x**2/sigma2) / np.sqrt(2*np.pi*sigma2)


#  ------------------------------------------  signal  -----------------------------------------  #


signal_freq = np.array((9.5,15.,18.,25.,31.))
signal_ampl = np.ones(5) * 0.2
signal_width = np.ones(5) * 0.1

def signalF(nu):
    temp = 0.
    
    for i in range(5):
        temp += signal_ampl[i] * gaussian(nu-signal_freq[i],signal_width[i]**2)
        
    return temp


nu = np.linspace(0.01,35,10000)
plt.plot(nu, signalF(nu), '-', label='signal')


#  ------------------------------------------  noise  ------------------------------------------  #

def noiseF(nu):
    # flat
    temp = 0.00119
    
    # gaussian of 13C
    omega = 2*np.pi*nu
    temp = 0.52 * np.exp(-0.5 * (omega-2.7156)*(omega-2.7156) * 1436)

    # 1/f
    if alpha>0:
        temp += A / nu**alpha
    
    return temp

nu = np.linspace(0.01,35,10000)
plt.plot(nu, noiseF(nu), '-', label='noise')


#  ----------------------------------------  spherical  ----------------------------------------  #
"""
filename = "Configurations/sSpher_T%.4f_dt%.4f_a%.4f.txt" % (Tfin,Delta_t,alpha)
filt = np.loadtxt(filename)

filtF = fft(filt, norm='forward')

nu = 0.5*np.pi*np.arange(N)/Tfin

plt.plot(nu[:N//2], np.abs(filtF[:N//2]), '-', label='spherical')
"""

#  ------------------------------------------  filter  -----------------------------------------  #

filename = "Configurations/s_T%.4f_dt%.4f_a%.4f_K%.4f_r%d.txt" % (Tfin,Delta_t,alpha,0,0)
filt = np.loadtxt(filename)

# find eta
f = open(filename, 'r')
head = f.readline()
f.close()
pos_equal = [i for i in range(len(head)) if head.startswith("=", i)]
etaInv = float( head[pos_equal[1]+1:-1] )

filtF = fftshift(fft(filt, norm='forward'))

nu = fftshift(fftfreq(N, d=Delta_t))
plt.plot(nu[N//2:], np.abs(filtF[N//2:]), '-', label='filter')
#plt.plot(nu[N//2:], np.real(filtF[N//2:]), '-', label='filter')

"""
ts = np.arange(N+1)*Delta_t
nu = np.linspace(0.01,35,1000)
filtF = np.zeros(len(nu))
for inu in range(len(nu)):
    #filtF[inu] = np.sum(filt*np.cos(2*np.pi*nu[inu]*np.arange(1,N+1)*Delta_t)) * Delta_t
    filtF[inu] = np.sum(filt * (np.sin(2*np.pi*nu[inu]*ts[1:]) - np.sin(2*np.pi*nu[inu]*ts[:-1])) ) / (2*np.pi*nu[inu])

plt.plot(nu, filtF, '-', label='filter')


filt = np.repeat(filt, 2)
filtF = fftshift(fft(filt, norm='forward'))

nu = fftshift(fftfreq(2*N, d=Delta_t/2))
plt.plot(nu[N:], np.abs(filtF[N:]), '-', label='filter')
"""
#  -------------------------------------------  plot  ------------------------------------------  #




plt.xlabel(r"$\nu$ [MHz]")
plt.ylabel(r"amplitude")

plt.xlim((-1,35))
plt.ylim((-0.01,0.9))

plt.title(r"$T$=%.2f $\mu$s, $\Delta t$=%.2f $\mu$s, $\alpha$=%.2f; $\eta^{-1}$=%.2f" % (Tfin,Delta_t,alpha, etaInv))

plt.legend(loc="upper right")
plt.show()

