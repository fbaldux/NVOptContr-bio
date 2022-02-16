/*  -------------------------------------------------------------------------------------------  //

    The program computes the field h for the spin glass Hamiltonian. The field represents
    the signal to be detected.
    The action potential of the neuron is modeled as a gaussian in time.

//  -------------------------------------------------------------------------------------------  */


#include <iostream>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;


// global variables
int N;
double Tfin, Delta_t;


//  -------------------------------  functions for the integral  ------------------------------  //

// integrated target signal
double integ_signal(double t_phase, double sigma, double t) {
    return sqrt(0.5*M_PI) * sigma * erf((t-t_phase)/(sqrt(2)*sigma));
    /*
    double x = (t-t_phase)/sigma;
    double temp = 1. + exp(x);
    temp = 0.5*exp(0.5*x)*(2.+temp) / (temp*temp) + 0.5*atan(exp(0.5*x));
    return temp * sqrt(sigma);
    */
}



//  ------------------------------------  save output data  -----------------------------------  //

void save_h(double *hs) {
    // create the output file
    char filename[100];
    snprintf(filename, 100, "Init/h_T%.4f_dt%.4f.txt", Tfin, Delta_t);        
    ofstream outfile(filename);
    
    if ( ! outfile.is_open() ) {
        cerr << "\nError with the output file!\n\n" << endl;
        exit(-1);
    }

    outfile << scientific;
    for (int k=0; k<N; k++) {
        outfile << hs[k] << endl;
    }
    
    outfile.close();
}


//  ------------------------------------------  main  -----------------------------------------  //

int main( int argc, char *argv[] ) {
    
    // parameter acquisition
    if( argc != 3 ) {
        cerr << "\nError! Usage: ./h_neuron <Tfin> <Delta_t>\n\n";
        exit(-1);
    }
    Tfin = strtof(argv[1], NULL);
    Delta_t = strtof(argv[2], NULL);

    // number of spins
    N = int(Tfin / Delta_t);
    
    
    // check if file already exists
    char filename[100];
    snprintf(filename, 100, "Init/h_T%.4f_dt%.4f.txt", Tfin, Delta_t); 
    ifstream outfile(filename);
    if (outfile) {
        exit(0);
    } 
    
    
    // dynamic allocation
    double *hs = new double[N]();   
    if (hs==NULL) {
        cerr << "\nError! Memory allocation failed.\n\n";
        exit(-1);
    }
    
    
    // integration: for the gaussian target it can be done analytically
    double sigma = 150;
    double t_phase = 400;
            
    for (int k=0; k<N; k++) {
        hs[k] = ( integ_signal(t_phase,sigma,(k+1)*Delta_t) - integ_signal(t_phase,sigma,k*Delta_t) ) / Tfin;
    }
    
    // save J to file
    save_h(hs);
    
    // free the memory
    free(hs);
    
    return 0;
}
