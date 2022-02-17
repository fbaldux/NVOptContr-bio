# OptControl-bio

Simulated annealing schedules for the optimal control problem of NV centers in diamond. Here, it is specified to target biological signals: the action potential of a neuron pulse, or the signal of a system of microtubules.



## anneal.sh

This shell runs the codes to

- compute the vectors `h[i]` and `J[i,j]` (for the experimental setup)
- perform many annealing cycles with `SA.cpp`, starting at high temperature from a random state. By default, it uses several values of K (the ferromagnetic coupling) parallelizing to different CPUs.


## exact+anneal.sh

This shell runs the codes to

- compute the vectors `h[i]` and `J[i,j]` (for the experimental setup)
- compute the exact solution with `spherical_FFT.py`
- perform some annealing cycles with `SA_from_spherical.cpp`, starting at low temperature from the exact solution.


## h\_microtub.cpp

The program computes the field h for the spin glass Hamiltonian. The field represents the signal to be detected.  
The microtubule field is represented by a 5-frequency signal.


## h\_neuron.cpp

The program computes the field h for the spin glass Hamiltonian. The field represents the signal to be detected.  
The action potential of the neuron is modeled as a gaussian in time.


## J\_C13\_1f.cpp

The program computes the couplings J for the spin glass Hamiltonian. The couplings represent the noise to be filtered out.  
The noise sources are the C13 impurities in the diamond, and a noise ~ amplNoise/f^alpha coming from the biological sample. Setting amplNoise=0 clearly removes the 1/f part of the noise.

NOTE: the C13 noise is

      S(omega) = S0 + A exp[ -(omega-omegaL)^2 / (2 sigma^2) ]

with

- S0 = 0.00119
- A = 0.52
- omegaL = 2 pi 0.4316 = 2.7118
- sigma  = 2 pi 0.0042 = 0.0264  =>  1/sigma^2 = 1434.8

NOTE: the 1/f noise is divergent at small frequencies, but the divergence is taken care of explicitly.


## J\_C13\_1f.py

Program that does the same of the .cpp version, just for checks.


## plot\_{...}.py

Just to plot the results.


## SA.cpp

The program anneals a random configuration of Ising spins s[i]=+/-1, according to the cost function

        H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - K sum_i s[i] s[i+1]

- The variables J[i,j] and h[i] are loaded from Init/
- The MC moves are only spin flips.
- The energy is computed efficiently at each step.
- The configurations found are saved to Configurations/, with the # of pulses and 1/eta in the header.
- The # of pulses and 1/eta for each configuration are saved to a file in Results/


## SA\_from\_spherical.cpp

The program anneals a random configuration of Ising spins s[i]=+/-1, according to the cost function

        H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - K sum_i s[i] s[i+1]

- The variables J[i,j] and h[i] are loaded from Init/
- The initial configuration is given by the output of spherical.py, i.e. the spherical model solution.
- The only allowed MC moves are domain wall shifts.
- The energy is computed efficiently at each step.
- The configurations found are saved to Configurations/, with the # of pulses and 1/eta in the header.
- The # of pulses and 1/eta for each configuration are saved to a file in Results/


## scatter.py

The program plots the sensitivities for all the values of K found in Results/


## spherical\_FFT.py

The program finds the configuration of _continuous_ spins `s[i]` that minimezes the cost function
   
      H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - lamda ( sum_i s[i]**2 - N )

- The variables `J[i,j]` and `h[i]` are loaded from Init/
- From the continuous spins are generated _Ising_ spins `s_Ising[i] = sign(s[i])`, that are then saved to Configurations/, with the # of pulses and 1/eta in the header.











