/*
 ===================================================================================================
 Name        : main.cpp
 Author      : Pasi Pyrrö
 Version     : 1.0
 Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis
 Date        : 21.6.2017
 Description : Sphere Decoder in C++14
 ===================================================================================================
 */

// #define ARMA_NO_DEBUG // for speed

#define _GLIBCXX_USE_CXX11_ABI 0

#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>

#include "misc.hpp"
#include "config.hpp"
#include "algorithms.hpp"
#include "sphdec.hpp"

using namespace std;
using namespace arma;

/* default filepaths */
string options_filename = "settings/settings.ini";
string basis_filename   = "bases/bases.txt";
string output_filename  = "output/output.txt";
string log_filename     = "logs/log.txt";

/* Storage for simulation parameters */
map<string, int> params;

/* random number generator */
mt19937_64 mersenne_twister{
    static_cast<long unsigned int>(
        chrono::high_resolution_clock::now().time_since_epoch().count()
    )
};

/* The program starts here */
int main(int argc, char** argv)
{

    string inputfile = options_filename;
    if (argc == 2){ // a parameter was given
        inputfile = "settings/" + string(argv[1]); // use alternative settings file
    } else if (argc > 2) { // too many parameters were given
        cout << "Usage: " << argv[0] << " [settings_file*]" << endl;
        exit(0);
    }

    configure(inputfile);

    log_msg("Starting program...");

    int m = params["no_of_transmit_antennas"];
    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];

    auto bases = read_matrices();
    cx_mat basis_sum(n, m);

    cout << "Read basis matrices: " << endl;
    for (auto const &basis : bases){
        cout << basis << endl;
        basis_sum += basis;
    }

    int *symbset = create_symbolset(params["x-PAM"]);
    
    log_msg("Using symbolset: " + vec2str(symbset, params["x-PAM"]));
    
    vector<pair<vector<int>,cx_mat>> codebook = create_codebook(bases, symbset);

    auto e = code_energy(codebook);
    log_msg("Code Energy");
    log_msg("-----------");
    log_msg("Average: " + to_string(e.first));
    log_msg("Max: " + to_string(e.second));

    

    double Hvar = 1, Nvar = 1;  

    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];

    

    /* initialize a bunch of complex matrices used in the simulation */
    cx_mat H, HX, X, N, Y, Ynorm;

    mat B(2*t*n, k), Q, R;

    vec y;
    vector<int> x(k);

    int runs = 0;
    int max_runs = 10000;
    int a = 0;
    int errors = 0; 
    // int max_err = params["required_errors"];

    double sigpow = 0;
    double noisepow = 0;
    double SNRreal = 0;
    double C = 0.0; // initial squared radius for the sphere decoder

    uniform_int_distribution<int> random_code(0, codebook.size()-1);

    log_msg("-----------");
    log_msg("Starting simulation...");

    // log_msg("Simulating received code matrices...");

    /* simulation main loop */
    for (int snr = min; snr <= max; snr += step) {
        // Hvar = e.first/pow(10, snr/10)*t; 
        Hvar = pow(10.0, snr/10.0)*(t/e.first); // calculate noise variance from SNR
        // cout << Hvar << endl; 
        while (/*errors < max_err && */ runs < max_runs){
            a = random_code(mersenne_twister);
            X = codebook[a].second; // Code block we want to send
            H = create_random_matrix(n, m, 0, Hvar);     // Channel matrix
            N = create_random_matrix(n, t, 0, Nvar);     // Noise matrix 

            // cout << X << endl << endl;
            // cout << H << endl << endl;
            // cout << N << endl << endl;
            // cout << "-------------" << endl;

            HX = H*X;
            sigpow += frob_norm_squared(HX);
            noisepow += frob_norm_squared(N);
            C = noisepow + 1e-3; // initial radius for the sphere decoder (added small "epsilon" to avoid equality comparison)
            // log_msg("Signal power: " + to_string(sigpow) + ", Noise power: " + to_string(noisepow));
            Y = HX + N; // Simulated code block that we would receive from MIMO-channel
            Ynorm = (Y + H*basis_sum*(q-1))*0.5; // normalize received matrix for the sphere decoder
            y = to_real_vector(Ynorm); // convert Y to real vector

            // B = (HX1 HX2 ... HXk)
            for(int i = 0; i < k; i++){
                B.col(i) = to_real_vector(H*bases[i]);
            }

            qr_econ(Q, R, B); // QR-decomposition of B (omits zero rows in R)
            y = Q.st()*y; // Map y to same basis as R

            x = sphdec(C, y, R, bases); // sphere decoder algorithm

            for (int j = 0; j < k; j++)     
                x[j] = 2*x[j] - q + 1;

            if (codebook[a].first != x){
                errors++;
            }

            // log_msg("Found point: " + vec2str(x, k) + ", sent point: " + vec2str(codebook[a].first, k));


            // errors += 1000;
            runs++;
        }
        
        // errors = 0;
        SNRreal = 10 * log(sigpow / noisepow) / log(10.0);
        noisepow = 0;
        sigpow = 0;
        
        log_msg("Simulated SNR: " + to_string(snr) + \
                ", Real SNR: " + to_string(SNRreal) + \
                ", BLER: " + to_string(errors) + "/" + to_string(runs) + " (" + to_string(100.0*(double)errors/runs) + " %)");

        runs = 0;
        errors = 0;
    }

    free(symbset);
    log_msg();
    return 0;
}