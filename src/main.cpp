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
#include "sphdec.hpp"
#include "algorithms.hpp"

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
    // log_msg("random error", "Error");
    // log_msg("random warning", "Warning");

    auto bases = read_matrices();
    log_msg("Read basis matrices:");
    for (auto const &base : bases){
        cout << base << endl;
    }

    int *symbset = create_symbolset(params["x-PAM"]);
    
    log_msg("Using symbolset: " + vec2str(symbset, params["x-PAM"]));
    
    auto codebook = create_codebook(bases, symbset);

    auto e = code_energy(codebook);
    log_msg("Code Energy");
    log_msg("-----------");
    log_msg("Average: " + to_string(e.first));
    log_msg("Max: " + to_string(e.second));

    cx_mat H, X, N, Y;

    double Hvar = 1, Nvar = 1;  

    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];

    int m = params["no_of_transmit_antennas"];
    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    // int k = params["no_of_matrices"];
    // int q = params["x-PAM"];


    int runs = 0;
    int max_runs = 1e5;
    // int errors = 0; 
    // int max_err = params["required_errors"];

    double sigpow = 0, noisepow = 0;
    double SNRreal = 0;

    // double C = 0.0; // initial squared radius for the sphere decoder

    uniform_int_distribution<int> random_code(0, codebook.size()-1);

    // log_msg("Simulating received code matrices...");

    /* simulation main loop */
    for (int snr = min; snr < max; snr += step) {
        Nvar = e.first/pow(10, snr/10); // calculate noise variance from SNR
        while (/*errors < max_err && */ runs < max_runs){
            X = codebook[random_code(mersenne_twister)]; // Code block we want to send
            H = create_random_matrix(n, m, 0, Hvar);     // Channel matrix
            N = create_random_matrix(n, t, 0, Nvar);     // Noise matrix 

            // C = frob_norm_squared(N) + 1e-3; // calculate initial radius

            sigpow += frob_norm_squared(H*X);
            noisepow += frob_norm_squared(N);

            Y = H*X + N; // Simulated code block that we would receive from MIMO-channel
            cout << Y << endl << endl;
            // errors += 1000;
            runs++;
        }
        // errors = 0;
        SNRreal = 10 * log(sigpow / noisepow) / log(10.0);
        noisepow = 0;
        sigpow = 0;
        log_msg("Simulated SNR: " + to_string(snr) + ", Real SNR: " + to_string(SNRreal));
    }

    free(symbset);
    log_msg();
    return 0;
}