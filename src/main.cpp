/*
 ===================================================================================================
 Name        : main.cpp
 Author      : Pasi Pyrrö
 Version     : 1.0
 Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis
 Date        : 18.7.2017
 Description : Sphere Decoder in C++14
 ===================================================================================================
 */

#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <csignal>
#include <thread>

#include "misc.hpp"
#include "config.hpp"
#include "algorithms.hpp"
#include "sphdec.hpp"

using namespace std;
using namespace arma;

/* storage for filepaths */
map<string, string> filenames;

/* Storage for simulation integer parameters */
map<string, int> params;

/* storage for simulation double parameters */
map<string, double> dparams;

/* storage for simulation string parameters */
map<string, string> sparams;

/* random number generator */
mt19937_64 mersenne_twister{
    static_cast<long unsigned int>(
        chrono::high_resolution_clock::now().time_since_epoch().count()
    )
};

/* Used for writing the csv output file */
parallel_vector<string> output;

/* flag which indicates early stop of the simulation */
bool exit_flag;

/* Handles task kills (CTRL-C) */
void signal_handler(int signum) {
    exit_flag = true;
    this_thread::sleep_for(chrono::milliseconds(1000));
    log_msg("Simulations terminated by user!", "Alert");
    output_data(output); // output current simulation data
    exit(signum);
}

/* The program starts here */
int main(int argc, char** argv)
{
    // assign signal SIGINT (when CTRL-C is pressed) to signal handler
    signal(SIGINT, signal_handler);
    exit_flag = false;

    /* define default filenames */
    filenames["settings"]         = "settings/settings.ini";
    filenames["settings_default"] = "settings/settings.ini";
    filenames["bases"]            = "bases/bases.txt";
    filenames["log"]              = "logs/log.txt";

    if (argc == 2){ // a parameter was given
        filenames["settings"] = "settings/" + string(argv[1]); // use alternative settings file
    } else if (argc > 2) { // too many parameters were given
        cout << "Usage: " << argv[0] << " [settings_file*]" << endl;
        exit(0);
    }

    configure();

    log_msg();
    log_msg("Starting program...");

    /* initialize simulation parameters */
    int m = params["no_of_transmit_antennas"];
    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];

    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];

    int min_runs = params["simulation_rounds"];

    int stat_interval = params["stat_display_interval"];

    double P = dparams["spherical_shaping_max_power"];

    string channel_model = sparams["channel_model"];

    vector<cx_mat> bases = read_matrices(filenames["bases"]);
    vector<cx_mat> coset_bases;

    bool coset_encoding = false;

    if (filenames.count("cosets") != 0)
        coset_encoding = true;

    if (coset_encoding)
        coset_bases = read_matrices(filenames["cosets"]);

    // cx_mat basis_sum(m, t);

    cout << "Read basis matrices: " << endl;
    for (auto const &basis : bases){
        cout << basis << endl;
        // basis_sum += basis;
    }
    // exit(0);

    cx_mat G = create_generator_matrix(bases);
    cx_mat invGe;

    if (coset_encoding)
        invGe = pinv(create_generator_matrix(coset_bases));

    cout << "Orthogonality check:" << endl;
    cout << G.t()*G << endl;

    mat Q, Rorig;
    qr_econ(Q, Rorig, create_real_generator_matrix(bases));  // QR-decomposition of G (omits zero rows in Rorig)
    process_qr(Q, Rorig);

    vector<int> symbset = create_symbolset(q);

    log_msg("Number of basis matrices (code length): " + to_string(k));
    log_msg("Using " + to_string(q) + "-PAM symbolset: " + vec2str(symbset, q));
    
    vector<pair<vector<int>,cx_mat>> codebook = create_codebook(bases, Rorig, symbset);

    auto e = code_energy(codebook);

    // e.first = 46.210487;
    
    log_msg("Code Energy");
    log_msg("-----------");
    log_msg("Average: " + to_string(e.first));
    log_msg("Max: " + to_string(e.second));
    log_msg("-----------");
    log_msg("Starting simulations... (Press CTRL-C to abort)");

    map<int, int> required_errors;
    if (filenames.count("error") == 0)
        for (int snr = min; snr <= max; snr += step)
            required_errors[snr] = params["required_errors"];
    else
        required_errors = read_error_requirements(filenames["error"]);

    if (filenames.count("output") == 0)
        create_output_filename();


    output.append("Simulated SNR,Real SNR,Runs,BLER,Avg Complexity"); // add label row
    
    #pragma omp parallel // parallelize SNR simulations
    {
        double Hvar = 1, Nvar = 1;
 
        /* initialize a bunch of complex matrices used in the simulation */
        cx_mat H, X(m, t), N(n, t); // Y //, HX, Ynorm;

        // mat B(2*t*n, k), Q, R, M;

        // vec y;
        // vec y2;

        vector<int> x(k), orig(k);
        // vec x(k);
        // vector<string> stats;

        int runs = 0;
        // int a = 0;
        int errors = 0; 
        int visited_nodes = 0;
        int total_nodes = 0;

        double sigpow = 0;
        double bler = 0;
        double avg_complex = 0;
        double noisepow = 0;
        double SNRreal = 0;
        double C = 0.0; // initial squared radius for the sphere decoder

        pair<vector<int>,cx_mat> codeword;

        // uniform_int_distribution<int> random_code(0, codebook.size()-1);

        /* simulation main loop */
        #pragma omp for
        for (int snr = min; snr <= max; snr += step) {
            // Hvar = e.first/pow(10, snr/10)*t; 
            Hvar = pow(10.0, snr/10.0)*(t/e.first); // calculate noise variance from SNR
            // cout << Hvar << endl; 
            while (errors < required_errors[snr] || runs < min_runs) {

                if (exit_flag) break; // terminate simulations

                // a = random_code(mersenne_twister);
                // X = codebook[a].second;                     
                if (P < 0)
                    codeword = create_random_codeword(bases, symbset); 
                else
                    codeword = create_random_spherical_codeword(bases, Rorig, symbset, P);

                orig = codeword.first;  // coefficients from the signal set (i.e. data vector)
                X = codeword.second;    // Code block we want to send

                // if (/*euclidean_norm(orig)*/ frob_norm_squared(X) > P + 1e-6 && P > 0) continue;  // We're outside the spherical constellation

                // cout << X << endl;
                // X = bases[0]*2.9+bases[1]*3;
                if (channel_model.compare("mimo") == 0)
                    H = create_random_matrix(n, m, 0, Hvar);    // Channel matrix
                else if (channel_model.compare("siso") == 0)
                    H = create_random_diag_matrix(n, 0, Hvar);
                else {
                    log_msg("Invalid channel model parameter used!", "Error");
                    exit(1);
                }
                // cout << H << endl;

                N = create_random_matrix(n, t, 0, Nvar);    // Noise matrix 
                // H.eye(2,2);
                // N.zeros(2,2);

                sigpow += frob_norm_squared(H*X);       // Signal power
                noisepow += frob_norm_squared(N);       // Noise power
                C = noisepow + 1e-3;                    // initial radius for the sphere decoder (added small "epsilon" to avoid equality comparison)

                // cout << N << endl;

                x = sphdec_wrapper(bases, Rorig, H, X, N, symbset, visited_nodes, C);

                // HX = H*X;
                // sigpow += frob_norm_squared(HX);
                // noisepow += frob_norm_squared(N);
                // C = noisepow + 1e-2; // initial radius for the sphere decoder (added small "epsilon" to avoid equality comparison)
                // // log_msg("Signal power: " + to_string(sigpow) + ", Noise power: " + to_string(noisepow));
                // Y = HX + N; // Simulated code block that we would receive from MIMO-channel
                // Ynorm = (Y + H*basis_sum*(q - 1))*0.5; // normalize received matrix for the sphere decoder
                // y = to_real_vector(Ynorm); // convert Y to real vector

                // // B = (HX1 HX2 ... HXk)
                // for(int i = 0; i < k; i++){
                //     B.col(i) = to_real_vector(H*bases[i]);
                // }
                // // cout << B << endl;

                // qr_econ(Q, R, B); // QR-decomposition of B (omits zero rows in R)
                // process_qr(Q, R); // Make sure R has positive diagonal elements
                // // cout << vec2str(y, y.n_elem) << endl;
                // y2 = Q.st()*y; // Map y to same basis as R


                // cout << "Input:" << endl << "C = " << C << endl;
                // cout << "y = " << vec2str(y2, y2.n_elem) << endl << endl;

                // cout << R << endl;
                // cout << Q << endl;
                // cout << "-----" << endl;

                // M.eye(4,4);
                // M *= 2;
                // // M(1,3) = 0.5;
                // // M(1,2) = 0.25;
                // // M(0,3) = -2.75;
                // // R(3,0) = 100;
                // // qr_econ(Q,R,R);
                // cout << "M = " << endl << M << endl;
                // // cout << Q << endl;
                // y2 = vec("0 0 0 0");
                // cout << "y = " << vec2str(y2, y2.size()) << endl;
                // for (int j = 0; j < k; j++)     
                //     y2[j] = 0.5*(y2[j] + q - 1);
                // cout << "y' = " << vec2str(y2, y2.size()) << endl;
                // y2 = M*y2;
                // cout << "M*y' = " << vec2str(y2, y2.size()) << endl;
                // C = 25;

                // #pragma omp task
                // x = R*Col<int>(sphdec(C, y2, R)); //, bases); // sphere decoder algorithm
                // x = sphdec(C, y2, R, visited_nodes);
                

                // if (x.size() == 0) { // point not found
                //     errors++;
                //     runs++;
                //     // Q.zeros();
                //     // R.zeros();
                //     continue;
                // }

                // for (int j = 0; j < k; j++)     
                //     x[j] = 2*x[j] - q + 1;
                    // orig[j] = 0.5*(codebook[a].first[j] + q - 1);

                // if (codebook[a].first != x){
                if(coset_encoding) {
                    if (!coset_check(G, invGe, Col<int>(orig) - Col<int>(x))) {
                        errors++;
                    }
                } else {
                    if (orig != x) {
                        errors++;
                    }
                }
                // cout << endl << vec2str(orig, orig.size()) << endl;
                
                // cout << "x = "<< vec2str(x, x.size()) << endl << endl;

                // log_msg("Found point: " + vec2str(x, k) + ", sent point: " + vec2str(codebook[a].first, k));

                total_nodes += visited_nodes;
                runs++;

                if (runs % stat_interval == 0 && stat_interval > 0){
                    SNRreal = 10 * log(sigpow / noisepow) / log(10.0);
                    bler = 100.0*(double)errors/runs;
                    avg_complex = (double)total_nodes/runs;
                    log_msg("SNR-simulation " + to_string(snr) + \
                    "\tReal SNR: " + to_string(SNRreal) + \
                    ", BLER: " + to_string(errors) + "/" + to_string(runs) + " (" + to_string(bler) + " %)" + \
                    ", Avg Complexity: " + to_string(avg_complex));
                    
                    // for (const string &s : stats)
                    //     log_msg(s);
                }

                // Q.zeros();
                // R.zeros();
                // x.clear();
                // orig.clear();
            }
            
            SNRreal = 10 * log(sigpow / noisepow) / log(10.0);
            bler = 100.0*(double)errors/runs;
            avg_complex = (double)total_nodes/runs;
            
            output.append(to_string(snr) + "," + to_string(SNRreal) + "," + to_string(runs) + "," + to_string(bler) + "," + to_string(avg_complex));
            
            log_msg("SNR-simulation " + to_string(snr) + \
                    "\t[Finished] Real SNR: " + to_string(SNRreal) + \
                    ", BLER: " + to_string(errors) + "/" + to_string(runs) + " (" + to_string(bler) + " %)" + \
                    ", Avg Complexity: " + to_string(avg_complex));

            // reset counters after simulation round
            runs = 0;
            errors = 0;
            noisepow = 0;
            sigpow = 0;
            total_nodes = 0;
            // stats.clear();
        }

    }
    /* output the simulation results in a csv file in /output/ folder */
    // if (output.size() > 1){
    //     sort(output.begin()+1, output.end(), snr_ordering); // sort vector by SNR (ascending order)
    //     log_msg("Printing simulation output to '" + filenames["output"] + "'...");
    //     output_csv(output);
    //     #ifdef PLOTTING // Draw plots with Gnuplot if plotting is enabled
    //     if (params["plot_results"] > 0){
    //         log_msg("Drawing plots...");
    //         plot_csv(1, 4, "SNR (dB)", "BLER (%)", true);
    //         plot_csv(1, 5, "SNR (dB)", "Average Complexity (# visited points)", false);
    //     }
    //     #endif
    // }
    output_data(output);

    log_msg("Program exited successfully!");
    // log_msg();
    return 0;
}