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

using namespace std;
using namespace arma;

string options_filename = "settings/settings.ini";
string basis_filename   = "bases/bases.txt";
string output_filename  = "output/output.txt";
string log_filename     = "logs/log.txt";

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
    cout << "Read basis matrices:" << endl;
    for (auto const &base : bases){
        cout << base << endl;
        // log_msg(base);
    }

    int *symbset = create_symbolset(params["x-PAM"]);
    cout << "Using symbolset: {";
    for (int i = 0; i < params["x-PAM"]; ++i) {
        cout << symbset[i];
        if (i < params["x-PAM"] - 1)
            cout << ", ";
    }
    cout << "}" << endl << endl;
    
    auto codebook = create_codebook(bases, symbset);

    auto e = code_energy(codebook);
    cout << "Code Energy" << \
    endl << "-----------" << \
    endl << "Average: " << e.first << \
    endl << "Max: " << e.second << endl;

    // cout << endl << "Random complex matrix test: " << endl << create_random_matrix(3,3,0,1) << endl;

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
    int errors = 0; 
    int max_err = params["required_errors"];

    uniform_int_distribution<int> random_code(0, codebook.size()-1);

    cout << endl << endl << "Simulating received code matrices..." << endl << endl;

    /* simulation main loop */
    for (int snr = min; snr < max; snr += step) {
        Nvar = e.first/pow(10, snr/10); // calculate noise variance from SNR
        while (errors < max_err){
            X = codebook[random_code(mersenne_twister)]; // Code block we want to send
            H = create_random_matrix(n, m, 0, Hvar);     // Channel matrix
            N = create_random_matrix(n, t, 0, Nvar);     // Noise matrix 

            Y = H*X + N; // Simulated received code block
            cout << Y << endl << endl;
            errors += 1000;
        }
        errors = 0;
    }

    free(symbset);
    log_msg();
    return 0;
}