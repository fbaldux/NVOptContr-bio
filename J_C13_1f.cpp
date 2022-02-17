/*  -------------------------------------------------------------------------------------------  //

    The program computes the couplings J for the spin glass Hamiltonian.
    The couplings represent the noise to be filtered out.

    The noise sources are the C13 impurities in the diamond, and a noise ~ amplNoise/f^alpha coming
    from the biological sample. Setting amplNoise=0 clearly removes the 1/f part of the noise. 

    NOTE: the C13 noise is
        S(omega) = S0 + A exp[ -(omega-omegaL)^2 / (2 sigma^2) ]
    with
        *) S0 = 0.00119
        *) A = 0.52
        *) omegaL = 2*pi*0.4316 = 2.7118
        *) sigma  = 2*pi*0.0042 = 0.0264  =>  1/sigma^2 = 1434.8

    NOTE: the 1/f noise is divergent at small frequencies, but the divergence is taken care of explicitly.

//  -------------------------------------------------------------------------------------------  */


#include <iostream>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

// global variables
int N;
double Tfin, Delta_t, alpha, amplNoise;
double invDt, invPi=1./M_PI;


//  -------------------------------  functions to be integrated  ------------------------------  //

// noise power spectrum: peaked contribution
double S1(double omega) {
    return 0.52 * exp(-0.5 * (omega-2.7118)*(omega-2.7118) * 1434.8);
}
// noise power spectrum: flat contribution
double S2(double omega) {
    return 0.00119;
}


// full integrand
// NOTE: x = omega*Delta_t
double integrand(double (*S)(double), double x, int k) {
    return S(x*invDt) * (1-cos(x)) * cos(k*x) / (x*x);
}

// integrand for the 1/f noise, with divergence at x=0 taken care of
// NOTE: omega = x/Delta_t
double integrand_1f(double x, int k) {
    return pow(x, -alpha) * ((1-cos(x)) * cos(k*x) / (x*x) - 0.5);
}


//  ----------------------------  numerical integration procedure  ----------------------------  //

// trapezoidal rule
double J_trapzInt(double (*S)(double), int k, double xMin, double xMax) {
    
    // integration steps: 1000 or at least 50 pts per cosine period
    int Nsteps = 50*k*(xMax-xMin) * (0.5*invPi);    
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
    
    return 4.*Delta_t*invPi * J*dx;
}


// integral of the function 1/x^alpha
double integrated_power(double xMin, double xMax) {
    if (alpha == 1.) {
        return log(xMax/xMin);
    }
    else {
        return (pow(xMax,1-alpha) - pow(xMin, 1-alpha)) / (1-alpha);
    } 
}


// trapezoidal rule handling the divergence at x=0 for the 1/f noise
double J_trapzInt_1f(int k, double xMin, double xMax) {
    
    // integration steps: 1000 or at least 50 pts per cosine period
    int Nsteps = 50*k*(xMax-xMin) * (0.5*invPi);    
    Nsteps = (Nsteps<1000) ? 1000 : Nsteps;
    
    double dx=(xMax-xMin)/Nsteps;
    
    // lower limit
    double J = 0.5*integrand_1f(xMin, k);
    // bulk
    for (int n=0; n<Nsteps; n++) {
        J += integrand_1f(xMin+n*dx, k);
    }
    // upper limit
    J += 0.5*integrand_1f(xMax, k);
    
    return 4.*amplNoise*invPi*pow(Delta_t,1.+alpha) * (J*dx + 0.5*integrated_power(xMin,xMax));
}


//  ------------------------------------  save output data  -----------------------------------  //

void save_J(double *Js) {
    // create the output file
    char filename[100];
    snprintf(filename, 100, "Init/J_T%.4f_dt%.4f_a%.4f_A%.2e.txt", Tfin, Delta_t, alpha, amplNoise); 
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
    if( argc != 5 ) {
        cerr << "\nError! Usage: ./J_1f <Tfin> <Delta_t> <alpha> <amplNoise>\n\n";
        exit(-1);
    }
    Tfin = strtof(argv[1], NULL);
    Delta_t = strtof(argv[2], NULL);
    invDt = 1./Delta_t;
    alpha = strtof(argv[3], NULL);
    amplNoise = strtof(argv[4], NULL);
    
    // number of spins
    N = int(Tfin / Delta_t);
    
    
    // check if file already exists
    char filename[100];
    snprintf(filename, 100, "Init/J_T%.4f_dt%.4f_a%.4f_A%.2e.txt", Tfin, Delta_t, alpha, amplNoise); 
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
        if (amplNoise > 0.) {
            Js[k] += J_trapzInt_1f(k, 0.001*Delta_t, 35*Delta_t);
        }
    }
    
    
    // save J to file
    save_J(Js);
    
    // free the memory
    free(Js);
    
    return 0;
}
