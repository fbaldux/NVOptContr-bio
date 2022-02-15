function mytime() {
    perl -MTime::HiRes=time -e 'printf "%.9f\n", time' 
    #echo "from time import time; print(time())"| python
}

Tfin=100.       # final time of the experiments
#Delta_ts=(0.5 0.2 0.1 0.05 0.02 0.01)       # pi-pulse distance
Delta_ts=( 1. )

annSteps=0       # number of steps in the temperature ramp
MCsteps=1       # number of MC steps at each ramp level
T0=0.01       # initial temperature
Reps=1       # number of states to sample

time0=$(mytime)

rm Init/* #Results/* Configurations/*

make J_experiment
make h_neuron
g++ -o SA SA_from_spherical.cpp -lm  -std=c++11

for i in {1..1}
do
    Delta_t=$Delta_ts[$i]
        
    ./J_experiment $Tfin $Delta_t
    #./h_neuron $Tfin $Delta_t 
    
    #python3 spherical_FFT.py $Tfin $Delta_t
    
    #./SA $Tfin $Delta_t $annSteps $MCsteps $T0 $Reps
    
done

rm J_experiment h_neuron SA





