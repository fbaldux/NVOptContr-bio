
rm Configurations/* Init/* Results/*

# time parameters
#Tfin=10.       # final time of the experiments
Delta_t=0.01       # pi-pulse distance

# noise parameters
#alpha=0

# annealing parameters
annSteps=1e3       # number of steps in the temperature ramp
MCsteps=1       # number of MC steps at each ramp level
T0=0.01       # initial temperature
Reps=1       # number of states to sample


g++ -o J J_C13_1f.cpp -lm
g++ -o h h_microtub.cpp -lm
g++ -o SA SA_from_spherical.cpp -lm -std=c++11


for Tfin in 20 #$(seq 1 1 16)
do
    ./h $Tfin $Delta_t
    
    for alpha in 1
    do
        #(
        ./J $Tfin $Delta_t $alpha
    
        python3 spherical_FFT.py $Tfin $Delta_t $alpha

        ./SA $Tfin $Delta_t $alpha $annSteps $MCsteps $T0 $Reps > E_${Tfin}_${Delta_t}_$alpha.txt
    
        python plot_Fourier.py $Tfin $Delta_t $alpha &
    
        #echo "Done Tfin=$Tfin"
        #)&
    done
    wait
done
rm J h SA

#python3 scatter.py $Tfin $Delta_t &




