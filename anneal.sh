
rm Configurations/*; rm Init/*; rm Results/*

# time parameters
Tfin=20.       # final time of the experiments
Delta_t=0.01       # pi-pulse distance

# noise parameters
alpha=1.5
amplNoise=1

# annealing parameters
annSteps=1e3       # number of steps in the temperature ramp
MCsteps=2e2       # number of MC steps at each ramp level
T0=0.1       # initial temperature
Reps=1       # number of states to sample

g++ -o J J_C13_1f.cpp -lm
g++ -o h h_microtub.cpp -lm
g++ -o SA SA.cpp -lm -std=c++11

./J $Tfin $Delta_t $alpha $amplNoise
#python J_C13_1f.py $Tfin $Delta_t $alpha $amplNoise &
./h $Tfin $Delta_t

for K in 1e-4 #5e-3 1e-3 5e-4 1e-4
do
    ./SA $Tfin $Delta_t $alpha $amplNoise $annSteps $MCsteps $T0 $K $Reps > Results/E_${Tfin}_${Delta_t}_$alpha.txt #&
    python verif.py $Tfin $Delta_t $alpha $amplNoise $K &
    :
done
wait 

rm J h SA


#python scatter.py $Tfin $Delta_t $alpha $amplNoise


