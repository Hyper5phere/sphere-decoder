/*
 ===================================================================================================
 Name        : main.cpp
 Author      : Pasi Pyrrö
 Version     : 1.0
 Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis
 Date        : 19.6.2017
 Description : Sphere Decoder in C++14
 ===================================================================================================
 */

// #define ARMA_NO_DEBUG // for speed

#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <regex>
#include <locale>
#include <random>
#include <chrono>
#include <mutex>
#include <ctime>

#include "misc.hpp"
#include "config.hpp"

using namespace std;
using namespace arma;

// string options_filename;
// string basis_filename;
// string output_filename;
// string log_filename;

string options_filename = "../settings/settings.ini";
string basis_filename   = "../bases/bases.txt";
string output_filename  = "../output/output.txt";
string log_filename     = "../logs/log.txt";

map<string, int> params;

// /* default filenames */
// string options_filename = "settings.ini";
// string basis_filename = "bases.txt";
// string output_filename = "output.txt";
// string log_filename = "log.txt";

// /* storage for simulation parameters */
// map<string, int> params;

/* random number generator */
mt19937_64 mersenne_twister{
    static_cast<long unsigned int>(
        chrono::high_resolution_clock::now().time_since_epoch().count()
    )
};

/* object for parallel computing synchronazation */
// mutex log_mutex;

/* function used for logging, should be thread safe */
// template <typename T>
// void log_msg(T msg){
//     struct tm *timeinfo;
//     time_t t = time(nullptr);
//     timeinfo = localtime(&t);
//     char timestamp[50]; 
//     strftime(timestamp, 50, "%d-%m-%Y %T | ", timeinfo);
//     cout << timestamp << msg << endl;
//     lock_guard<mutex> lock(log_mutex); // make sure other threads don't write to log file simultaneosly
//     ofstream logfile(log_filename, ios_base::app);
//     logfile << timestamp << msg << endl;
//     logfile.close();
// }

/* signum function */
int sign(double x){
    return (x < 0) ? -1 : 1;
}

/* Takes the squared Frobenius norm from a complex matrix A */
double frob_norm_squared(cx_mat A){
    double sum = 0;
    for (auto i = 0u; i < A.n_rows; i++)
        for (auto j = 0u; j < A.n_cols; j++)
            sum += norm(A(i,j)); 
    return sum;
}

/* generates a random n x m complex matrix from uniform distribution */
cx_mat create_random_matrix(int n, int m, double mean, double variance){
    uniform_real_distribution<double> distr(mean, variance);
    cx_mat A(n,m);
    return A.imbue([&]() {
        return complex<double>(distr(mersenne_twister), distr(mersenne_twister));
    });
}

/* Creates a C-style integer array representation of an x-PAM symbolset */
int* create_symbolset(int q = params["x-PAM"]){
    int *symbset = (int*) malloc(q * sizeof(int));
    for (int u = 0; u < q; u++) {
        symbset[u] = 2*u - q + 1;
    }
    return symbset; // must be free()'d after use
}

/* Calculates all combinations of elements for code vector _a_
   given the set of feasible symbols (element values) */
void combinations(set< vector<int> > &comblist, vector<int> symbset, vector<int> comb, int dim){
    comblist.insert(comb);
    if (dim >= 0){
        for (const int symbol : symbset){
            comb[dim] = symbol;
            combinations(comblist, symbset, comb, dim-1);
        }
    }
}

/* Helper function for above combinations algorithm */
set< vector<int> > comb_wrapper(int* symbset, int vector_len){
    set< vector<int> > comblist;
    vector<int> init(vector_len);
    vector<int> symbset_v(symbset, symbset + params["x-PAM"]);
    for (int i = 0; i < vector_len; i++)
        init[i] = symbset[0];
    combinations(comblist, symbset_v, init, vector_len-1);
    return comblist;
}

/* Creates a codebook (set of X matrices) from basis matrices B_i and symbolset x-PAM */
vector<cx_mat> create_codebook(const vector<cx_mat> &bases, int* symbolset){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    // int cs = params["codebook_size"];
    int samples = params["energy_estimation_samples"];

    /* lattice generator matrix G (alternative approach) */
    // cx_mat G(m*n,k);
    // for(int i = 0; i < k; i++){
    //     G.col(i) = vectorise(bases[i]);
    // }

    vector<cx_mat> codebook;
    cx_mat X(m, t);
    
    if(samples > 0){
        int random_index = 0;
        uniform_int_distribution<int> dist(0,q-1);
        cout << "Random sampled data vector combinations:" << endl;
        for (int j = 0; j < samples; j++){
            cout << "{";
            for (int i = 0; i < k; i++) {
                random_index = dist(mersenne_twister);
                cout << symbolset[random_index] << " ";
                X = X + symbolset[random_index]*bases[i];
            }
            cout << "}" << endl;
            codebook.push_back(X);
            X.zeros();
            // cout << endl << X << endl;

            // cout << "{";
            // for (int s = 0; s < k - 1; s++)
            //     cout << symbolset[s] << ", ";
            // cout << symbolset[k-1] << "}" << endl;
        }
        cout << endl;
    } else {

        /* all possible combinations of code words */
        auto c = comb_wrapper(symbolset, k);

        cout << "All possible data vector combinations:" << endl;
        for (const auto &symbols : c){
            for (int i = 0; i < k; i++)
                X = X + symbols[i]*bases[i];
            
            codebook.push_back(X);
            X.zeros();

            cout << "{";
            for (int s = 0; s < k - 1; s++)
                cout << symbols[s] << ", ";
            cout << symbols[k-1] << "}" << endl;
        }
        cout << endl;
    }
    return codebook;
}

/* Computes the average and maximum energy of given codebook X */
pair<double,double> code_energy(const vector<cx_mat> X){
    double sum = 0, max = 0, tmp = 0, average = 0;
    int cs = (int) X.size();

    for (int i = 0; i < cs; i++){
        //tmp = pow(norm(X[i], "fro"), 2);
        tmp = frob_norm_squared(X[i]);
        sum += tmp;
        // cout << X[i] << endl;
        // cout << tmp << endl;
        // cout << sum << endl;
        if (tmp > max)
            max = tmp;
    }
    average = sum / cs;
    return make_pair(average, max);
}

/* The program starts here */
int main(int argc, char** argv)
{
    // options_filename = "../settings/settings.ini";
    // basis_filename   = "../bases/bases.txt";
    // output_filename  = "../output/output.txt";
    // log_filename     = "../logs/log.txt";

    string inputfile = options_filename;
    if (argc == 2){ // a parameter was given
        inputfile = "../settings/" + string(argv[1]); // use alternative settings file
    } else if (argc > 2) { // too many parameters were given
        cout << "Usage: " << argv[0] << " [settings_file*]" << endl;
        exit(0);
    }

    configure(inputfile);

    // misc::log_msg("[Info] Program has started.");

    auto bases = read_matrices();
    cout << "Read basis matrices:" << endl;
    for (auto const &base : bases){
        cout << base << endl;
        // log_msg(base);
    }

    int *symbset = create_symbolset();
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
    // cout << endl << "sign() test: " << endl << sign(-100.0) << sign(19) << sign(0) << endl;

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


    int errs = 0, errors = params["required_errors"];

    uniform_int_distribution<int> random_code(0, codebook.size()-1);

    cout << endl << endl << "Simulating received code matrices..." << endl << endl;

    /* simulation main loop */
    for (int snr = min; snr < max; snr += step) {
        Nvar = e.first/pow(10, snr/10); // calculate noise variance from SNR
        while (errs < errors){
            X = codebook[random_code(mersenne_twister)]; // Code block we want to send
            H = create_random_matrix(n, m, 0, Hvar);     // Channel matrix
            N = create_random_matrix(n, t, 0, Nvar);     // Noise matrix 

            Y = H*X + N; // Simulated received code block
            cout << Y << endl << endl;
            errs += 1000;
        }
        errs = 0;
    }

    free(symbset);
    return 0;
}