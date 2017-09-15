/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : main.cpp                                                                              *
 * Project     : Schnorr-Euchnerr sphere decoder simulation for space-time lattice codes               *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++11                                                                                 *
 * Description : main program file, handles the simulation setup and main loop                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 
#define ARMA_NO_DEBUG /* disable Armadillo bound checks for addiotional speed */

#include <iostream>
#include <iomanip>
#include <armadillo> /* linear algebra library */
#include <complex>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <chrono>
#include <ratio>
#include <ctime>
#include <csignal>
#include <thread>

#include "misc.hpp"
#include "config.hpp"
#include "algorithms.hpp"
#include "sphdec.hpp"

using namespace std;
using namespace arma;

/* storage for filepaths */
unordered_map<string, string> filenames;

/* Storage for simulation integer parameters */
unordered_map<string, int> params;

/* storage for simulation double parameters */
unordered_map<string, double> dparams;

/* storage for simulation string parameters */
unordered_map<string, string> sparams;

/* 64-bit standard random number generator */
mt19937_64 mersenne_twister {
    static_cast<long unsigned int>(
        chrono::high_resolution_clock::now().time_since_epoch().count() /* seed for rng */
    )
};

/* general purpose large uniform integer distribution (use with modulo) */
uniform_int_distribution<int> uniform_dist(0, 1000000);
uniform_real_distribution<double> uniform_real_dist(0.0, 1.0);
// normal_distribution uniform_dist(0, 1000000);

/* Used for writing the csv output file */
parallel_vector<string> output;

/* flag which indicates early stop of the simulation */
bool exit_flag;

/* Handles task kills (CTRL-C) */
void signal_handler(int signum) {
    exit_flag = true;    /* Terminates simulations */
    this_thread::sleep_for(chrono::milliseconds(1000)); /* Give some time for the simulation to end 
                                                           before attempting to output data */
    log_msg("Simulations terminated by user!", "Alert");
    output_data(output); /* output current simulation data */
    exit(signum);        /* exit program */
}

/* The program starts here */
int main(int argc, char** argv)
{
    /* assign signal SIGINT (when CTRL-C is pressed) to signal handler */
    signal(SIGINT, signal_handler);
    /* if this flag is true, the program will exit from simulation loops */
    exit_flag = false;

    /* define default filenames */
    filenames["settings"]         = "settings/settings.ini";
    filenames["settings_default"] = "settings/settings.ini";
    filenames["bases"]            = "bases/bases.txt";
    filenames["log"]              = "logs/log.txt";

    /* settings filename parameter was given */
    if (argc == 2){ 
        filenames["settings"] = "settings/" + string(argv[1]); 
    }
    /* too many parameters were given */
    else if (argc > 2) { 
        cout << "Usage: " << argv[0] << " [settings_file*]" << endl;
        exit(0);
    }

    /* This program uses custom logging function which prints info
       in stdout (console output) and also in a log file.
       It can be found from the misc.cpp file */
    log_msg();
    log_msg("Starting program...");

    /* Configure program parameters i.e. occupy the maps listed above 
     * with values read from the settings file.
     * Implemented in config.cpp file
     */
    configure();

    /* initialize simulation parameters */
    int m = params["no_of_transmit_antennas"]; 
    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int s = params["codebook_size_exponent"]; /* if spesified, the program will 
                                                 try to look for a spherical shaping radius 
                                                 so that the codebook has size of 2^s */
    
    /* Spesify SNR range for the simulation */
    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];

    int min_runs = params["simulation_rounds"]; /* run at least this many rounds for each SNR simulation */

    int stat_interval = params["stat_display_interval"];  /* Display intermediate simulation stats every this many rounds */
    int search_density = params["radius_search_density"]; /* Spesifies how many radiuses around the initial guess we try
                                                             to get the codebook size to 2^s */ 
    search_density = (search_density <= 0) ? 10 : search_density; /* default to 10 */

    double P = dparams["spherical_shaping_max_power"]; /* User can also manually spesify the spherical shaping radius */

    string channel_model = sparams["channel_model"]; /* Decides how channel matrix H is generated */

    /* Read the complex basis matrices for the lattice we want to decode in */
    vector<cx_mat> bases = read_matrices(filenames["bases"]);
    /* Another set of basis matrices in case we want to use coset encoding */
    vector<cx_mat> coset_bases;

    /* flag that spesifies if we're doing the so called 'wiretap' simulation 
       i.e. all codewords that belong to the same coset have the same meaning 
       for this simulation we need another lattice, which is a sublattice of the previous one */
    bool coset_encoding = false;

    /* If the user has spesified a coset basis file, assume 'wiretap' simulation */
    if (filenames.count("cosets") != 0)
        coset_encoding = true;

    /* Conditionally read the coset basis matrices from a different basis file */
    if (coset_encoding)
        coset_bases = read_matrices(filenames["cosets"]);

    /* Sanity check: Print all read matrices! */
    log_msg("Read basis matrices:");
    for (auto const &basis : bases){
        log_msg(mat2str(basis), "Raw");
    }
    if (coset_encoding) {
        log_msg("Read coset basis matrices:");
        for (auto const &basis : coset_bases){
            log_msg(mat2str(basis), "Raw");
        }
    }

    /* Generator matrix of the lattice code */
    cx_mat G = create_generator_matrix(bases);

    /* LLL_reduction for G */
    // cout << G << endl;
    // G = LLL_reduction(G);
    // cout << G << endl;
    // bases = generator_to_bases(G);

    /* output LLL reduced basis matrices */
    // auto lll_bases = generator_to_bases(G);
    // for (const auto &lll : lll_bases) {
    //     output_complex_matrix("bases/MIDO_LLL_basis.txt", lll, true);
    // }
    // exit(0);

    /* sublattice generator matrices for coset encoding (e stands for 'Eve') */
    cx_mat invGe, Ge;

    /***** If the Gram matrix (generator times its own transpose) is diagonal, 
           the basis for the lattice is orthogonal, Check it like this *****/
    // cout << "Orthogonality check:" << endl;
    // cout << G.t()*G << endl;

    /***** Matrices can be hardcoded in the simulation like so *****/
    // mat coset_multiplier("4 0 0 0;"
    //                      "0 2 0 0;"
    //                      "0 0 2 0;"
    //                      "0 0 0 2");

    // mat coset_multiplier("-2 -2  0  0;"
    //                       "0  0 -2 -1;"
    //                      "-1  1  1 -2;"
    //                       "1 -1  1 -1");

    // mat coset_multiplier("4 2 2 2;"
    //                      "0 2 0 0;"
    //                      "0 0 2 0;"
    //                      "0 0 0 2");

    // mat coset_multiplier("4 0 0 0;"
    //                      "0 4 0 0;"
    //                      "0 0 4 0;"
    //                      "0 0 0 4");

    // mat coset_multiplier("-2 -3  4 -1;"
    //                      " 0 -1  0  3;"
    //                      " 0 -3 -2 -3;"
    //                      "-4 -1  0 -1");

    if (coset_encoding) {
        /* as q-PAM is already a subset of translated 2Z^n lattice, 
           we need to multiply the sublattice basis with 2 to make the lattice points comparable */
        Ge = 2*create_generator_matrix(coset_bases); /* Read cosets from file */
        // Ge = 2*G*coset_multiplier; /* Generate cosets from basis matrix with hardcoded transformation matrices */
        invGe = pinv(Ge); /* Coset lattice generator matrix and its (pseudo)inverse */
    }

    /***** You can generate basis matrix txt files like this *****/
    // auto cos_bases = generator_to_bases(0.5*Ge);
    // for (const auto &b : cos_bases) {
    //     output_complex_matrix("bases/alamouti_L5.txt", b, true);
    // }

    /***** Another example of matrix outputting: prints the real generator matrix into a txt file *****/
    // output_real_matrix("bases/MIDO_basis.txt", G_real);

    /* Many helper functions have been implemented for matrix manipulation,
       Check algorithms.cpp for a complete list */
    mat G_real = to_real_matrix(G); /* Make a real representation out of G */

    /***** Calculate the lattice constant like this *****/    
    // cout << "determinant: " << det(G_real.t()*G_real) << endl;

    mat Q, Rorig; /* Decompose G_real */
    qr_econ(Q, Rorig, G_real);  /* QR-decomposition of G_real (omits zero rows in Rorig) */
    process_qr(Q, Rorig); /* Makes sure Rorig only has positive diagonal elements (this is probably optional) */

    /* Creates a q-PAM symbol set,
       e.g. 4-PAM = {-3, -1, 1, 3} 
       Details in algorithms.cpp */
    vector<int> symbset = create_symbolset(q);

    // cout << "floor 6.5: " << nearest_symbol_floor(6.5, symbset) << endl;
    // cout << "floor -10: " << nearest_symbol_floor(-10, symbset) << endl;
    // cout << "floor  10: " << nearest_symbol_floor(10, symbset) << endl;
    // cout << "ceil  5.5: " << nearest_symbol_ceil(5.5, symbset) << endl;
    // cout << "ceil  -10: " << nearest_symbol_ceil(-10, symbset) << endl;
    // cout << "ceil   10: " << nearest_symbol_ceil(10, symbset) << endl;


    int num_points = 0; /* Number of codewords inside the spherical constellation */
    
    /* try to find radius for the codebook that has atleast 2^s codewords if s is spesified */
    if (s > 0) {
        log_msg("Attempting to estimate squared radius (max power) for codebook of 2^" + to_string(s) + " codewords...");
        double P_estimate = estimate_squared_radius(Rorig, s); /* see algorithms.cpp for details */
        cx_vec lambda_min = shortest_basis_vector(G); /* pick the shortest column vector of G */
        double search_step = frob_norm_squared(lambda_min)/search_density; /* radius search step */
        double rcurr = P_estimate + search_density*search_step/2; /* starting radius */
        log_msg("Initial guesstimate for codebook squared radius: " + to_string(P_estimate));
        log_msg("Radius search step: " + to_string(search_step));
        vector<double> rvec; /* store radius candidates in here (descending order) */
        for (int i = 0; i < search_density + 1 && rcurr > 0.0; i++) {
            rvec.push_back(rcurr);
            rcurr -= search_step;
        }
        /* Counts corresponding number of points inside hypersphere for each radius in rvec */
        vector<int> pvec = count_points_many_radiuses(Rorig, symbset, rvec, vec(k, fill::zeros), k, 0);
        /* Pick the smallest radius that maps to number of codewords in range [2^s, 2^(s+1)] */
        for (auto j = pvec.size()-1; j >= 0; j--) {
            if (pvec[j] >= (int)pow(2, s) && pvec[j] < (int)pow(2, s+1)) {
                P = rvec[j];
                cout << P << endl;
                num_points = pvec[j];
                break;
            }
        }
        /* Evaluate estimation result */
        if (P < 0) {
            log_msg("Estimation failed, using non-spherical shaping...", "Alert");
        } else if (P == dparams["spherical_shaping_max_power"]) {
            log_msg("Estimation failed, using the configured value for squared radius...", "Alert");
        } else {
            dparams["spherical_shaping_max_power"] = P;
        }
    }

    /* Generate the whole codebook or a random sampled subset of it. */
    log_msg("Generating codebook...");
    vector<pair<vector<int>,cx_mat>> codebook = create_codebook(bases, Rorig, symbset);

    /* Calculate the average and maximum energy of the codewords in our lattice code */
    auto e = code_energy(codebook);

    /* Print some useful information before starting the actual simulation */
    log_msg("", "Raw");
    log_msg("Simulation info");
    log_msg("---------------");
    log_msg("Number of basis matrices (code length): " + to_string(k));
    log_msg("Using " + to_string(q) + "-PAM symbolset: " + vec2str(symbset, q));
    log_msg("Codebook size: " + to_string(codebook.size()));
    log_msg("Average code energy: " + to_string(e.first));
    log_msg("Max code energy: " + to_string(e.second));
    if (P > 0) {
        log_msg("Using codebook spherical shaping squared radius: " + to_string(P));
        // log_msg("Suggested squared radius (max power) for 2^" + to_string(s) +
        //     " (" + to_string((int)pow(2, s)) + ") codewords: " + to_string(P_estimate));
        if (num_points == 0)
            num_points = count_points(Rorig, symbset, P, vec(k, fill::zeros), k, 0);
        log_msg("Number of codewords inside the hypersphere: " + to_string(num_points));
    }
    if (coset_encoding) {
        auto rates = code_rates(2*G, Ge);  /* Calculate code rates, scaling by two needed again because of q-PAM (Ge already scaled) */
        log_msg("Code overall rate: "      + to_string(get<0>(rates)/t) + " bpcu");
        log_msg("Code transmission rate: " + to_string(get<1>(rates)/t) + " bpcu");
        log_msg("Code confusion rate: "    + to_string(get<2>(rates)/t) + " bpcu");
    }
    log_msg("---------------");
    log_msg("", "Raw");

    /* Ask the user whether he/she actually wants to run the simulations after precalculations */
    string answer;
    cout << "Continue to simulation (y/n)? " << endl;
    getline(cin, answer);
    /* Continue only if user types 'y' or 'yes' */
    if (!(answer.compare("y") == 0 || answer.compare("yes") == 0)) {
        log_msg("Program exited successfully!");
        return 0;
    }

    log_msg("Starting simulations... (Press CTRL-C to abort)");

    unordered_map<int, int> required_errors;
    /* if error file is spesified, fill the required_errors map with it
       Otherwise use constant error requirement for each SNR simulation */
    if (filenames.count("error") == 0)
        for (int snr = min; snr <= max; snr += step)
            required_errors[snr] = params["required_errors"];
    else
        required_errors = read_error_requirements(filenames["error"]);

    int snr_to_errorbound = 0;

    /* If user has not spesified output filename, generate one automatically */
    if (filenames.count("output") == 0)
        create_output_filename(); /* see misc.cpp for details */

    /* thread safe output string vector, used for csv output */
    output.append("Simulated SNR,Real SNR,Avg Complexity,Max Complexity,Errors,Runs,BLER"); /* add label row */

    log_msg("------------------------------------------------------------------------------------------------------------");
    log_msg("SNR-simulation | Real SNR     | Avg Complexity | Max Complexity | Errors     | Runs       | BLER ");
    log_msg("------------------------------------------------------------------------------------------------------------");

    auto start = chrono::high_resolution_clock::now();
    
    #pragma omp parallel /* parallelize SNR simulations */
    {
        double Hvar = 1, Nvar = 1;
 
        /* initialize a bunch of complex matrices used in the simulation */
        cx_mat H(n, m), X(m, t), N(n, t);

        vector<int> x(k), orig(k); /* output and input vectors */

        ostringstream buffer; /* string stream buffer used for log_msg() */

        /* simulation variables */
        int runs = 0;
        int errors = 0; 
        int visited_nodes = 0;
        int total_nodes = 0;
        int samples = params["energy_estimation_samples"]; /* Take this many random samples from the codebook (negative number means sample all) */

        double sigpow = 0;
        double bler = 0;
        double avg_complex = 0;
        int    max_complex = 0;
        double noisepow = 0;
        double SNRreal = 0;
        double C = 0.0; /* initial squared radius for the sphere decoder */

        /* pair containing the coefficients and matrix representation of a single codeword */
        pair<vector<int>,cx_mat> codeword;
        vector<pair<vector<int>,cx_mat>> codebook_instance = codebook;
        uniform_int_distribution<int> dist(0, codebook.size()-1);

        /* Simulations main loop (iteration shared among parallel threads) 
           Simulations are indexed by their SNR value */
        #pragma omp for schedule(static,1)
        for (int snr = min; snr <= max; snr += step) {

            Hvar = pow(10.0, snr/10.0)*(t/e.first); /* calculate channel matrix variance from SNR */
            /* (noise variance is constant one) */

            /* Single SNR simulation loop: 
               run until both conditions are satisfied or the simulation is terminated by user */
            snr_to_errorbound = required_errors[snr];
            while (errors < snr_to_errorbound || runs < min_runs) {

                if (exit_flag) break; /* terminate simulations */
                 
                /* Simulation starts by generating a random codeword that we want to send */
                if (samples <= 0) {
                    codeword = codebook_instance[dist(mersenne_twister)];
                } else {
                    if (P <= 0){
                        codeword = create_random_codeword(bases, symbset); 
                    }
                    else {
                        codeword = create_random_spherical_codeword(bases, Rorig, symbset, P);
                    }
                }


                orig = codeword.first;  /* coefficients from the signal set (i.e. data vector) */
                X = codeword.second;    /* Code block we want to send */

                /* Generate random channel matrix according to channel mode */
                if (channel_model.compare("mimo") == 0)
                    H = create_random_matrix(n, m, 0, Hvar);
                else if (channel_model.compare("siso") == 0)
                    H = create_random_diag_matrix(n, 0, Hvar);
                else {
                    log_msg("Invalid channel model parameter used!", "Error");
                    exit(1);
                }

                N = create_random_matrix(n, t, 0, Nvar); /* Additive complex Gaussian white noise matrix */

                sigpow += frob_norm_squared(H*X);       /* Signal power */
                noisepow += frob_norm_squared(N);       /* Noise power */
                C = noisepow + 1e-3;                    /* initial radius for the sphere decoder (added small "epsilon" to avoid equality comparison) */

                /* wrapper function for the sphere decoder algorithm, see details in sphdec.cpp 
                   x is the decoded output vector */
                x = sphdec_wrapper(bases, Rorig, H, X, N, symbset, visited_nodes, C);

                /* Check if the decoded vector is the same as what we sent in the simulation. 
                   If we're doing 'wiretap' simulation the checking is done in a different manner */
                if (coset_encoding) {
                    /* see details for coset_check in algorithms.cpp 
                       the last argument is the difference vector between x and orig */
                    if (!coset_check(G, invGe, Col<int>(orig) - Col<int>(x))) {
                        errors++;
                    }
                } else {
                    /* simple equality check for normal simulation */
                    if (orig != x) {
                        errors++;
                    }
                }

                /***** the vectors can be printed out nicely like this for debugging reasons *****/
                // cout << endl << "orig = " << vec2str(orig, orig.size()) << endl;
                // cout << "x = "<< vec2str(x, x.size()) << endl << endl;

                /* increase counters */
                total_nodes += visited_nodes;
                if (visited_nodes > max_complex) max_complex = visited_nodes;
                runs++;

                /* print intermediate simulation stats every now and then (if enabled) */
                if (runs % stat_interval == 0 && stat_interval > 0){
                    SNRreal = 10 * log(sigpow / noisepow) / log(10.0);
                    bler = (double)errors/runs;
                    avg_complex = (double)total_nodes/runs;
                    /*log_msg("SNR-simulation " + to_string(snr) + \
                    "\tReal SNR: " + to_string(SNRreal) + \
                    ", BLER: " + to_string(errors) + "/" + to_string(runs) + " (" + to_string(bler) + ")" + \
                    ", Avg Complexity: " + to_string(avg_complex) + \
                    ", Max Complexity: " + to_string(max_complex));*/

                    buffer << std::right << std::setw(14) << snr << " | " << setw(12) << SNRreal << " | "
                           << std::setw(14) << avg_complex << " | " << std::setw(14) << max_complex << " | " 
                           << std::setw(10) << errors << " | " << std::setw(10)  << runs << " | " << std::setw(16) << bler;

                    log_msg(buffer.str());
                    buffer.str("");
                }
            }
            
            /* SNR simulation finished, calculate results... */

            /* Sanity check for the realized SNR, should roughly equal to the configured SNR */
            SNRreal = 10 * log(sigpow / noisepow) / log(10.0); 
            /* Block error rate */
            bler = (double)errors/runs;
            /* Average complexity (number of visited search tree nodes) of the sphere decoder */
            avg_complex = (double)total_nodes/runs;
            
            output.append(to_string(snr) + "," + float2str(SNRreal, 16) + "," + float2str(avg_complex, 16) + "," + to_string(max_complex)
                  + "," + to_string(errors) + "," + to_string(runs) + "," + float2str(bler, 16));
            

            /*log_msg("SNR-simulation " + to_string(snr) + \
                    "\t[Finished] Real SNR: " + to_string(SNRreal) + \
                    ", BLER: " + to_string(errors) + "/" + to_string(runs) + " (" + to_string(bler) + ")" + \
                    ", Avg Complexity: " + to_string(avg_complex) + \
                    ", Max Complexity: " + to_string(max_complex));*/


            buffer << "[Finished] " << std::right << std::setw(3) << snr << " | " << setw(12) << SNRreal << " | "
                   << std::setw(14) << avg_complex << " | " << std::setw(14) << max_complex << " | " 
                   << std::setw(10) << errors << " | " << std::setw(10)  << runs << " | " << std::setw(16) << bler;

            log_msg(buffer.str());
            buffer.str("");

            /* reset counters after simulation round */
            runs = 0;
            errors = 0;
            noisepow = 0;
            sigpow = 0;
            total_nodes = 0;
        }

    }
    log_msg("------------------------------------------------------------------------------------------------------------");
    auto end = chrono::high_resolution_clock::now();
    auto simulation_time = chrono::duration_cast<chrono::duration<double>>(end - start);
    log_msg("Simulations finished in " + to_string(simulation_time.count()) + " seconds.");

    /* output the simulation results in a csv file in /output/ folder */
    output_data(output);

    log_msg("Program exited successfully!");
    return 0;
}