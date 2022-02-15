function mytime() {
    perl -MTime::HiRes=time -e 'printf "%.9f\n", time' 
    #echo "from time import time; print(time())" | python
}

#Tfin=4.       # final time of the experiments
Delta_t=0.02       # pi-pulse distance

annSteps=1e3       # number of steps in the temperature ramp
MCsteps=1       # number of MC steps at each ramp level
T0=0.01       # initial temperature
Reps=1       # number of states to sample

time0=$(mytime)


g++ -o J J_experiment.cpp -lm
g++ -o h h_microtub.cpp -lm
g++ -o SA SA_from_spherical.cpp -lm -std=c++11

for Tfin in $(seq 10 10 300) #1 2 3 4 5 10 20 30 40 50 100
do
    ./J $Tfin $Delta_t
    ./h $Tfin $Delta_t
    
    python3 spherical_FFT.py $Tfin $Delta_t #&& echo "spherical done" $( echo "$(mytime) - $time0" | bc -l )

    ./SA $Tfin $Delta_t $annSteps $MCsteps $T0 $Reps
done
rm J h SA

#python3 scatter.py $Tfin $Delta_t &




