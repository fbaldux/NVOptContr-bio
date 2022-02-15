# OptControl-SimAnneal

Simulated annealing schedule for the optimal control problem of NV centers in diamond (here specified to target the action potential of a neuron pulse). 


---
---
## Programs

---
### h\_neuron.cpp

The program computes the field h for the spin glass Hamiltonian. The field represents the signal to be detected.  
Currently, the supported option is a gaussian in time.


---
### J\_experiment.cpp

The program computes the couplings J for the spin glass Hamiltonian. The couplings represent the noise to be filtered out.


---
### SA.cpp

The program anneals a random configuration of Ising spins s[i]=+/-1, according to the cost function

        H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - K sum_i s[i] s[i+1]

- The variables J[i,j] and h[i] are loaded from Init/
- The MC moves are only spin flips.
- The energy is computed efficiently at each step.
- The configurations found are saved to Configurations/, with the # of pulses and 1/eta in the header.
- The # of pulses and 1/eta for each configuration are saved to a file in Results/


---
### SA\_from\_spherical.cpp

The program anneals a random configuration of Ising spins s[i]=+/-1, according to the cost function

        H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - K sum_i s[i] s[i+1]

- The variables J[i,j] and h[i] are loaded from Init/
- The initial configuration is given by the output of spherical.py, i.e. the spherical model solution.
- The only allowed MC moves are domain wall shifts.
- The energy is computed efficiently at each step.
- The configurations found are saved to Configurations/, with the # of pulses and 1/eta in the header.
- The # of pulses and 1/eta for each configuration are saved to a file in Results/


---
### spherical\_FFT.py

The program finds the configuration of _continuous_ spins `s[i]` that minimezes the cost function
   
      H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - lamda ( sum_i s[i]**2 - N )

- The variables `J[i,j]` and `h[i]` are loaded from Init/
- From the continuous spins are generated _Ising_ spins `s_Ising[i] = sign(s[i])`, that are then saved to Configurations/, with the # of pulses and 1/eta in the header.


---
### spherical\_diag.py

The program finds the configuration of _continuous_ spins `s[i]` that minimezes the cost function
   
      H = 0.5 sum_ij J[i,j] s[i] s[j] - log |sum_i h[i] s[i]| - lamda ( sum_i s[i]**2 - N )

- The variables `J[i,j]` and `h[i]` are loaded from Init/
- Exact diagonalization is used instead of the FFT.
- From the continuous spins are generated _Ising_ spins `s_Ising[i] = sign(s[i])`, that are then saved to Configurations/, with the # of pulses and 1/eta in the header.


---
---
## Shells

---
### anneal.sh

This shell runs the codes to

- compute the vectors `h[i]` and `J[i,j]` (for the experimental setup)
- perform many annealing cycles with `SA.cpp`, starting at high temperature from a random state. By default, it uses several values of K (the ferromagnetic coupling) parallelizing to different CPUs.


---
### exact.sh

This shell runs the codes to

- compute the vectors `h[i]` and `J[i,j]` (for the experimental setup)
- compute the exact solution with `spherical_FFT.py`


---
### exact+anneal.sh

This shell runs the codes to

- compute the vectors `h[i]` and `J[i,j]` (for the experimental setup)
- compute the exact solution with `spherical_FFT.py`
- perform some annealing cycles with `SA_from_spherical.cpp`, starting at low temperature from the exact solution.











