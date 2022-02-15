function mytime() {
    perl -MTime::HiRes=time -e 'printf "%.9f\n", time' 
    #echo "from time import time; print(time())" | python
}

rm Configurations/* Init/* Results/*

# time parameters
Tfin=16.       # final time of the experiments
Delta_t=0.02       # pi-pulse distance

# noise parameters
alpha=1.5

# annealing parameters
annSteps=1e3       # number of steps in the temperature ramp
MCsteps=1       # number of MC steps at each ramp level
T0=0.01       # initial temperature
Reps=1       # number of states to sample

time0=$(mytime)

#if g++ -o J J_experiment.cpp -lm
if g++ -o J J_1f.cpp -lm
then
    ./J $Tfin $Delta_t $alpha
    echo "J done" $( echo "$(mytime) - $time0" | bc -l )
fi


#if g++ -o h h_neuron.cpp -lm
if g++ -o h h_microtub.cpp -lm
then
    ./h $Tfin $Delta_t
    echo "h done" $( echo "$(mytime) - $time0" | bc -l )
fi


python3 spherical_FFT.py $Tfin $Delta_t $alpha && echo "spherical done" $( echo "$(mytime) - $time0" | bc -l )


if g++ -o SA SA_from_spherical.cpp -lm -std=c++11
then
    ./SA $Tfin $Delta_t $alpha $annSteps $MCsteps $T0 $Reps &
    wait
    
    echo "SA done" $( echo "$(mytime) - $time0" | bc -l )
fi
rm J h SA

#python3 scatter.py $Tfin $Delta_t &

python noisePlot.py $Tfin $Delta_t $alpha &



