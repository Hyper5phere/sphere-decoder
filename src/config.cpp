/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : config.cpp                                                                            *
 * Project     : Planewalker - Schnorr-Euchner sphere decoder simulation for space-time lattice codes  *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++ (2011 or newer standard)                                                          *
 * Description : Implements functions for handling user input and simulation configuration             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ARMA_NO_DEBUG /* disable Armadillo bound checks for addiotional speed */

#include <iostream>
#include <armadillo> /* linear algebra library */
#include <complex>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>

#include "config.hpp"
#include "misc.hpp"

#define NUM_OPTIONS 22

using namespace std;
using namespace arma;

/* Creates a default configuration file in case one does not exist */
void create_config(){
    ofstream defconf(filenames["settings_default"]);
    // the file can contain comments similar to this line here
    defconf << "// configuration settings and simulation parameters for the sphere decoder program" << endl << endl;
    defconf << "basis_file=bases.txt            // Text file containing the basis matrices" << endl;
    defconf << "coset_file=                     // Optionally specify the coset encoding sublattice basis matrix file" << endl;
    defconf << "output_file=                    // Optionally specify the output filename" << endl;
    defconf << "error_file=                     // Optionally specify the file containing error requirements for the SNR simulations" << endl;
    defconf << "channel_model=mimo              // Define the channel model for the simulation (either 'mimo' or 'siso')" << endl;
    defconf << "x-PAM=4                         // The size of the PAM signaling set (even positive integer)" << endl;
    defconf << "energy_estimation_samples=-1    // Number of samples to make the code energy estimation (-1 = sample all)" << endl;
    defconf << "no_of_matrices=2                // Number of basis matrices (dimension of the data vectors)" << endl;
    defconf << "matrix_coefficient=1.0          // Multiply all basis matrices by this constant" << endl;
    defconf << "time_slots=2                    // Number of time slots used in the code" << endl;
    defconf << "no_of_transmit_antennas=2       // Number of transmit antennas" << endl;
    defconf << "no_of_receiver_antennas=2       // Number of receiver antennas"  << endl;
    defconf << "snr_min=6                       // Minimum value for signal-to-noise ratio" << endl;
    defconf << "snr_max=12                      // Maximum value for signal-to-noise ratio" << endl;
    defconf << "snr_step=2                      // Increase SNR by this value per each iteration" << endl;
    defconf << "simulation_rounds=100000        // Number of simulation rounds to run" << endl;
    defconf << "required_errors=500             // Demand at minimum this many errors before the simulation ends" << endl;
    defconf << "plot_results=-1                 // Draw plots? (1 = yes, -1 = no)" << endl;
    defconf << "stat_display_interval=-1        // Defines after each how many rounds to display the current simulation stats (-1 = disabled)" << endl;
    defconf << "spherical_shaping_max_power=-1  // Defines the maximum distance from origin for codebook elements (-1 = unbounded)" << endl;
    defconf << "codebook_size_exponent=-1       // The codebook will have 2^s codewords where s is this parameter (overrides above parameter)" << endl;
    defconf << "radius_search_density=-1        // Defines how accurate the codebook squared radius estimation will be" << endl;
    defconf.close(); 
}

/* Read simulation parameters from a settings file */
void configure() {
    string filepath = filenames["settings"];
    string default_filepath = filenames["settings_default"];

    ifstream config_file(filepath);

    if (!config_file.good() && filepath.compare(default_filepath) != 0) {
        log_msg("No settings file '" + filepath + "' found, using the default one...", "Alert");
        config_file = ifstream(default_filepath);
    } 

    if (!config_file.good()) {
        log_msg("No default settings file found, creating a new one with default settings...");
        create_config();
        log_msg("Make your changes to the settings file and rerun this program.");
        log_msg("Exiting...");
        exit(0);
    }

    /* List special (not integer) parameter names that should be handled differently */
    vector<string> sparam_names = {"channel_model"};
    vector<string> dparam_names = {"spherical_shaping_max_power", "matrix_coefficient"};
    vector<string> accept_empty_names = {"output_file", "error_file", "coset_file"};
    
    string line;
    char *endptr; /* used for parsing error checking */
    int lines = 0;
    while (getline(config_file, line))
    {
        istringstream iss(line);
        string key;
        if(getline(iss, key, '=') )
        {
            if (key.find("//", 0) != string::npos)
                continue;
            lines++;
            // cout << lines << endl;
            string value;
            getline(iss, value);
            clean_input(value); 
            if (!value.empty()) { /* Valid key-value pair found */
                /* kinda messy checking but works */
                if (key.compare("basis_file") == 0) { 
                    filenames["bases"] = "bases/" + value;
                } else if (key.compare("coset_file") == 0) {
                    filenames["cosets"] = "bases/" + value;
                } else if (key.compare("output_file") == 0) {
                    filenames["output"] = "output/" + value;
                } else if (key.compare("error_file") == 0) {
                    filenames["error"] = "settings/" + value;
                } else if (find(dparam_names.begin(), dparam_names.end(), key) != dparam_names.end()) {
                    dparams[key] = strtod(value.c_str(), &endptr);
                    if (*endptr != 0) {
                        log_msg("Invalid value for option '" + key + "'", "Error");
                        exit(1);
                    }
                } else if (find(sparam_names.begin(), sparam_names.end(), key) != sparam_names.end()) {
                    sparams[key] = value;
                } else {
                    params[key] = strtol(value.c_str(), &endptr, 10); 
                    if (*endptr != 0) {
                        log_msg("Invalid value for option '" + key + "'", "Error");
                        exit(1);
                    }
                }
            } else if (find(accept_empty_names.begin(), accept_empty_names.end(), key) == accept_empty_names.end()){
                log_msg("[Error] Value for option '" + key + "' not spesified!");
                exit(1);
            }
        }
    }
    if (lines < NUM_OPTIONS){
        log_msg("Too few options spesified in the '" + filepath + "'!", "Error");
        log_msg("Consider deleting the default settings file which will reset the program settings.");
        exit(1);
    } 
    if ((params["x-PAM"] <= 0) || (params["x-PAM"] % 2 == 1)){
        log_msg("Invalid x-PAM signal set spesified!", "Error");
        log_msg("x-PAM option accepts positive even integers.");
        exit(1);
    }

    // helper variable (size of the codebook) calculated from the input parameters
    dparams["codebook_size"] = pow(params["x-PAM"], params["no_of_matrices"]);
}

/* reads k (m x t) matrices from the spesified basis_file */
vector<cx_mat> read_matrices(const string &filepath){

    /* - uses regular expression pattern to search 
     *   matrix elements in almost any known format (eg. Mathematica, Matlab (with commas))
     * - ignores most white spaces and supports asterisk in front of 'I'
     * - pattern is also case insensitive
     * - Caveats: elements can't be separated by spaces alone like in Matlab, 
     *            only one sign (x/-) is to be inserted between real and imag part of a complex number
     */

    /* Some supported formats for the matrix elements:
     * 1) "re ,"
     * 2) "re + im *I,"
     * 3) "re - im *I,"
     * 4) "im *I," 
     */

    /* The regex pattern */
    regex elem_regex(
        "("                                                                                      // three basic types for matrix element
            "([+-]?\\s*[^\\d\\.]+[Ii]{1})"                                                       // 1. plain letter 'i' with sign
            "|" 
            "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)\\s*\\*?[Ii]?)"                       // 2. pure real or complex number
            "|"
             "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)"                                     // 3. both real part
            "\\s*"                                                                               
            "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*|[+-]?\\s*[^\\d\\.]+)\\s*\\*?[Ii]{1})" // and imaginary part
        ")" 
        "\\s*[,)};\\]\n]{1}"       // element ends in comma, semicolon, newline or closing bracket (WHITE SPACE SEPARATION NOT SUPPORTED)
    );

    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int idx = 0;

    ifstream matrix_file(filepath);
    vector<cx_mat> output;
    smatch match; /* regex match object, contains capture groups */
    vector< complex<double> > numbers; /* store the parsed elements in here */
    
    string content((istreambuf_iterator<char>(matrix_file)), 
        istreambuf_iterator<char>());

    /* capture group substrings */
    string  m2,  // plain i
            m3,  // real or pure complex number
            m5,  // real part
            m6,  // imaginary part with i
            m7;  // imaginary part without i
    
    /* Search the configured basis_file for basis matrix elements */
    while(regex_search(content, match, elem_regex)){

        /***** print all capture groups for debugging! *****/
        // for (auto i = 0u; i < match.size(); ++i)
        //     cout << i << ": " << match[i].str() << endl;
        // cout << endl;

        /* assign and clean the helper substrings */
        m2 = match[2].str();
        m3 = match[3].str();
        m5 = match[5].str();
        m6 = match[6].str();
        m7 = match[7].str();

        clean_input(m2);
        clean_input(m3);
        clean_input(m5);
        clean_input(m6);
        clean_input(m7);

        /* both real and imaginary part are present */
        if (!m5.empty() && !m6.empty()) {
            if (isdigit(m7.back()) || (m7.back() == '.')) 
                numbers.push_back(complex<double>(stod(m5), stod(m7)));
            else { /* imag element is of type +-i */
                numbers.push_back(complex<double>(stod(m5), stod(m7 + string("1.0"))));
            }
        } 
        /* plain +-i element found */
        else if (!m2.empty()) {
            remove_from_string(m2, 'i');
            remove_from_string(m2, 'I');
            numbers.push_back(complex<double>(0.0, stod(m2 + string("1.0"))));
        }
        /* pure real or imag element found */
        else {
            if (isdigit(m3.back()) || (m3.back() == '.')) {
                numbers.push_back(complex<double>(stod(m3), 0.0)); 
            } else {
                numbers.push_back(complex<double>(0.0, stod(m3)));
            }
        }
        /* iterate the content, i.e. remove the part of the string we examined already */
        content = match.suffix().str();
        idx++;
    }

    /* Sanity check */
    if (m*t*k != (int)numbers.size()){
        log_msg("Failed to read the " + to_string(m*t*k) + " matrix elements (found " + to_string(numbers.size()) + \
            ") configured in '" + filepath + "'!", "Error");
        exit(1);
    }

    /* store the read matrix elements in a complex matrix datatype from Armadillo library */
    idx = 0;
    cx_mat X(m, t);
    for (int i = 0; i < k; i++){    
        for (int j = 0; j < m; j++){
            for (int s = 0; s < t; s++){
                X(j,s) = numbers[idx];
                idx++;
            }
        }
        output.push_back(dparams["matrix_coefficient"]*X);
        X.zeros();
    }
    
    return output;
}

/* Reads the user spesified error requirements for each SNR simulation from a csv file 
 * Example: 'settings/errors.csv' file contains two lines:
 * --------------------
 * 10,12,14,16,18
 * 1000,800,600,400,200
 * --------------------
 * first line has SNR values separated by commas and 
 * second line has the corresponding error requirement for each SNR value
 * this means that SNR 16 simulation will be run until we hit 400 errors etc.
 */
unordered_map<int, int> read_error_requirements(const string &filepath){
    unordered_map<int, int> snr_to_error;
    vector<int> snrs, errors;
    int line_num = 1;
    ifstream error_file(filepath);
    string line;
    string value;

    if (!error_file.good()) {
        log_msg("No error requirements file '" + filepath + "' found!", "Error");
        exit(1);
    }

    while(getline(error_file, line)){
        clean_input(line);
        if (line.empty()) continue;
        istringstream iss(line);
        while(getline(iss, value, ',')){
            // clean_input(value);
            if (line_num == 1){
                snrs.push_back(strtol(value.c_str(), NULL, 10));
            } else if (line_num == 2){
                errors.push_back(strtol(value.c_str(), NULL, 10));
            } else break;
        }
        line_num++;
    }

    /* assert error requirements validity... */
    /* especially make sure that the data here matches with the data in the settings file */

    if (snrs.empty()){
        log_msg("Error requirements file '" + filepath + "' is empty!", "Error");
        exit(1);
    }

    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];
    unsigned required_size = (max-min)/step + 1;

    if (snrs.size() != required_size){
        log_msg("Error requirements file '" + filepath + "' has an incorrect amount of SNR values!", "Error");
        exit(1);
    }

    if (snrs.size() != errors.size()){
        log_msg("Error requirements file '" + filepath + "' has an incorrect amount of error values!", "Error");
        exit(1);
    }

    for (auto i = 0u; i < snrs.size(); i++){
        if (min + (int)i*step != snrs[i]){
            log_msg("Error requirements file '" + filepath + "' has unspesified SNR value at column " + to_string(i+1) + "!", "Error");
            exit(1);
        }
        snr_to_error[snrs[i]] = errors[i];
    }

    /* All good! */

    return snr_to_error;
}