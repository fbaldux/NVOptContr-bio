function mytime() {
    perl -MTime::HiRes=time -e 'printf "%.9f\n", time' 
    #echo "from time import time; print(time())" | python
}

Tfin=100.       # final time of the experiments
Delta_t=0.02       # pi-pulse distance

annSteps=5e3       # number of steps in the temperature ramp
MCsteps=2e2       # number of MC steps at each ramp level
T0=0.1       # initial temperature
Reps=5       # number of states to sample

time0=$(mytime)

if g++ -o J J_experiment.cpp -lm
then
    ./J $Tfin $Delta_t
    echo "J done" $( echo "$(mytime) - $time0" | bc -l )
fi


#if g++ -o h h_neuron.cpp -lm
if g++ -o h h_microtub.cpp -lm
then
    ./h $Tfin $Delta_t
    echo "h done" $( echo "$(mytime) - $time0" | bc -l )
fi


if g++ -o SA SA.cpp -lm -std=c++11
then
    for K in 5e-3 1e-3 5e-4 1e-4
    do
        ./SA $Tfin $Delta_t $annSteps $MCsteps $T0 $K $Reps &
    done
    wait 
    echo "SA done" $( echo "$(mytime) - $time0" | bc -l )
fi
rm J h SA


python scatter.py $Tfin $Delta_t


