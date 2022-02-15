/*  -------------------------------------------------------------------------------------------  //

    The program computes the field h for the spin glass Hamiltonian. The field represents
    the signal to be detected.
    
    Currently, the supported option is a gaussian in time.

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

// integrated monochromatic wave
double integ_mono(double A, double nu0, double t) {
    return 0.5/M_PI * A * sin(2*nu0*M_PI*t) / nu0;
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
    
    
    // integration: equal weight superposition of certain frequencies
    double A[5]  = {0.2,0.2,0.2,0.2,0.2};
    double nu[5] = {9.5,15.,18.,25.,31.};
    
    for (int k=0; k<N; k++) {
        // it should be already initialized but you never know
        hs[k] = 0.;
        
        for (int j=0; j<5; j++) {
            hs[k] += ( integ_mono(A[j],nu[j],(k+1)*Delta_t) - integ_mono(A[j],nu[j],k*Delta_t) ) / Tfin;
        }
    }
    
    // save J to file
    save_h(hs);
    
    // free the memory
    free(hs);
    
    return 0;
}
