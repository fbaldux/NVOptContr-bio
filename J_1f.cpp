/*  -------------------------------------------------------------------------------------------  //

    The program computes the couplings J for the spin glass Hamiltonian. The couplings 
    represent the noise to be filtered out.
    
//  -------------------------------------------------------------------------------------------  */


#include <iostream>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

// global variables
#define amplNoise 0.5
int N;
double Tfin, Delta_t, invDt, alpha;


//  -------------------------------  functions to be integrated  ------------------------------  //

// noise power spectrum: peaked contribution
double S1(double omega) {
    return 0.52 * exp(-0.5 * (omega-2.7156)*(omega-2.7156) * 1436);
}
// noise power spectrum: flat contribution
double S2(double omega) {
    return 0.00119;
}
// noise power spectrum: 1/f^alpha contribution
double S3(double omega) {
    return amplNoise * 2*M_PI / omega;
}


// full integrand: x = omega*Delta_t
double integrand(double (*S)(double), double x, int k) {
    return S(x*invDt) * (1-cos(x)) * cos(k*x) / (x*x);
}


//  ----------------------------  numerical integration procedure  ----------------------------  //

// trapezoidal rule
double J_trapzInt(double (*S)(double), int k, double xMin, double xMax) {
    
    // integration steps: 500 or at least 50 pts per cosine period
    int Nsteps = 50*k*(xMax-xMin) / (2*M_PI);
    Nsteps = (Nsteps<1000) ? 1000 : Nsteps;
        
    double dx=(xMax-xMin)/Nsteps;
    
    // lower limit
    double J = 0.5*integrand(S, xMin, k);
    // bulk
    for (int n=0; n<Nsteps; n++) {
        J += integrand(S, xMin+n*dx, k);
    }
    // upper limit
    J += 0.5*integrand(S, xMax, k);
    
    return 4*Delta_t/M_PI * J*dx;
}



//  ------------------------------------  save output data  -----------------------------------  //

void save_J(double *Js) {
    // create the output file
    char filename[100];
    snprintf(filename, 100, "Init/J_T%.4f_dt%.4f_a%.4f.txt", Tfin, Delta_t, alpha); 
    ofstream outfile(filename);
    
    if ( ! outfile.is_open() ) {
        cerr << "\nError with the output file!\n\n" << endl;
        exit(-1);
    }

    outfile << scientific;
    for (int k=0; k<N; k++) {
        outfile << Js[k] << endl;
    }
    
    outfile.close();
}


//  ------------------------------------------  main  -----------------------------------------  //

int main( int argc, char *argv[] ) {
    
    // parameter acquisition
    if( argc != 4 ) {
        cerr << "\nError! Usage: ./J_1f <Tfin> <Delta_t> <alpha>\n\n";
        exit(-1);
    }
    Tfin = strtof(argv[1], NULL);
    Delta_t = strtof(argv[2], NULL);
    invDt = 1./Delta_t;
    alpha = strtof(argv[3], NULL);
    
    // number of spins
    N = int(Tfin / Delta_t);
    
    
    // check if file already exists
    char filename[100];
    snprintf(filename, 100, "Init/J_T%.4f_dt%.4f_a%.4f.txt", Tfin, Delta_t, alpha); 
    ifstream outfile(filename);
    if (outfile) {
        exit(0);
    } 
    
    
    // dynamic allocation
    double *Js = new double[N]();   
    if (Js==NULL) {
        cerr << "\nError! Memory allocation failed.\n\n";
        exit(-1);
    }
    
    
    // numerical integration
    for (int k=0; k<N; k++) {
        Js[k] = J_trapzInt(S1, k, 2.4517*Delta_t, 2.9795*Delta_t);
        Js[k] += J_trapzInt(S2, k, 0.001*Delta_t, 6.5*Delta_t);
        if (alpha > 0.) {
            Js[k] += J_trapzInt(S3, k, 0.001*Delta_t, 10*Delta_t);
        }
    }
    
    
    // save J to file
    save_J(Js);
    
    
    // free the memory
    free(Js);
    
    return 0;
}
