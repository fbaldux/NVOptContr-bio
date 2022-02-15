Delta_t=0.16
tone=3       # monochromatic or trichromatic
harmonic=0       # number of the harmonic (only for the monochromatic signal)

annSteps=1e3       # number of steps in the temperature ramp
Reps=10       # number of states to sample


make J_experiment
make h_experiment
g++ -o SA SA.cpp -lm -std=c++11
g++ -o SA_from_spherical SA_from_spherical.cpp -lm  -std=c++11

echo "# Tfin t_anneal t_exact t_ann_from_exact"  > times.txt

for (( Tfin=50.; Tfin<=400.; Tfin+=25. )) 
do  
    echo -n "$Tfin " >> times.txt
    
    ./J_experiment $Tfin $Delta_t
    ./h_experiment $Tfin $Delta_t $tone $harmonic
    
    ./SA $Tfin $Delta_t $tone $harmonic $annSteps 1e2 0.1 0.003 $Reps 
    
    python3 spherical_FFT.py $Tfin $Delta_t $tone $harmonic >> Results/cont_T$Tfin.txt
    ./SA_from_spherical $Tfin $Delta_t $tone $harmonic $annSteps 1 0.01 $Reps 
    
    echo '' >> times.txt
    echo "$Tfin done"
done
rm J_experiment h_experiment SA SA_from_spherical





