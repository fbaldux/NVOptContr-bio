
#rm Configurations/*; rm Init/*; rm Results/*

# time parameters
#Tfin=10.       # final time of the experiments
Delta_t=0.01       # pi-pulse distance

# noise parameters
#alpha=0
amplNoise=1

# annealing parameters
annSteps=1e3       # number of steps in the temperature ramp
MCsteps=1       # number of MC steps at each ramp level
T0=0.01       # initial temperature
Reps=1       # number of states to sample


g++ -o J J_C13_1f.cpp -lm
g++ -o h h_microtub.cpp -lm
g++ -o SA SA_from_spherical.cpp -lm -std=c++11


for Tfin in $(seq 2 2 16)
do
    ./h $Tfin $Delta_t
    
    for alpha in 0.5 1 1.5
    do
        (
        ./J $Tfin $Delta_t $alpha $amplNoise
        #python J_C13_1f.py $Tfin $Delta_t $alpha $amplNoise 
    
        python3 spherical_FFT.py $Tfin $Delta_t $alpha $amplNoise

        ./SA $Tfin $Delta_t $alpha $amplNoise $annSteps $MCsteps $T0 $Reps #> Results/E_${Tfin}_${Delta_t}_$alpha.txt
    
        #python plot_Fourier.py $Tfin $Delta_t $alpha $amplNoise 0 &    
        )&
    done
    wait
    
    echo "Done Tfin=$Tfin"
done
rm h SA

#python3 scatter.py $Tfin $Delta_t &




