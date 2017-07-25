#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <regex>

#include "config.hpp"
#include "misc.hpp"

#define NUM_OPTIONS 18

using namespace std;
using namespace arma;

/* Creates a default configuration file in case one does not exist */
void create_config(){
    ofstream defconf(filenames["settings_default"]);
    // the file can contain comments similar to this line here
    defconf << "// configuration settings and simulation parameters for the sphere decoder program" << endl << endl;
    defconf << "basis_file=bases.txt            // Text file containing the basis matrices" << endl;
    defconf << "output_file=auto                // Optionally spesify the output filename (auto = automatic)" << endl;
    defconf << "channel_model=mimo              // define the channel model for the simulation (either 'mimo' or 'siso')" << endl;
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
        // log_msg();
        exit(0);
    }

    vector<string> sparam_names = {"channel_model"};
    vector<string> dparam_names = {"spherical_shaping_max_power", "matrix_coefficient"};
    
    string line;
    int lines = 0;
    while (getline(config_file, line))
    {
        istringstream iss(line);
        string key;

        if(getline(iss, key, '=') )
        {
            if (key.find("//", 0) != string::npos)
                continue;
            string value;
            getline(iss, value);
            clean_input(value); 
            if (!value.empty()) {
                // cout << value << endl;
                if (key.compare("basis_file") == 0) {
                    filenames["bases"] = "bases/" + value;
                } else if (key.compare("output_file") == 0) {
                    if (value.compare("auto") != 0)
                        filenames["output"] = "output/" + value;
                } else if (find(dparam_names.begin(), dparam_names.end(), key) != dparam_names.end()) {
                    if ((dparams[key] = strtod(value.c_str(), NULL)) == 0) {
                        log_msg("Invalid value for option '" + key + "'", "Error");
                        // log_msg();
                        exit(1);
                    }
                } else if (find(sparam_names.begin(), sparam_names.end(), key) != sparam_names.end()) {
                    sparams[key] = value;
                } else {
                    if ((params[key] = strtol(value.c_str(), NULL, 10)) == 0) {
                        log_msg("Invalid value for option '" + key + "'", "Error");
                        // log_msg();
                        exit(1);
                    }
                }
            } else {
                log_msg("[Error] Value for option '" + key + "' not spesified!");
                // log_msg();
                exit(1);
            }
        }
        lines++;
    }
    if (lines < NUM_OPTIONS){
        log_msg("Too few options spesified in the '" + filepath + "'!", "Error");
        log_msg("Consider deleting the default settings file which will reset the program settings.");
        // log_msg();
        exit(1);
    } 
    if ((params["x-PAM"] <= 0) || (params["x-PAM"] % 2 == 1)){
        log_msg("Invalid x-PAM signal set spesified!", "Error");
        log_msg("x-PAM option accepts positive even integers.");
        // log_msg();
        exit(1);
    }

    // helper variable (size of the code matrix set) calculated from the input parameters
    // params["codebook_size"] = (int)pow(params["x-PAM"], params["no_of_matrices"]);
}

/* reads k (m x t) matrices from the spesified basis_file */
vector<cx_mat> read_matrices(){

    /* - regular expression pattern used to search 
     *   matrix elements in almost any known format (eg. Mathematica, Matlab)
     * - ignores most white spaces and supports asterisk in front of 'I'
     * - pattern is also case insensitive
     */
    regex elem_regex(
        //"(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*\\*?[Ii]?\\s*[,)};\\]\n]{1})|"
        "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*\\s*\\*?[Ii]?)|"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*" // white space
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*\\*?[Ii]{1})\\s*([,)};\\]\n]{1})" // element ends in comma, semicolon, newline or closing bracket
    );

    /* Some supported formats for the matrix elements:
     * 1) "re ,"
     * 2) "re + im *I,"
     * 3) "re - im *I,"
     * 4) "im *I," 
     */

    // Used to determine whether the element has both real and complex part
    regex complex_split("\\d+\\.?\\s*[+-]{1}");

    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int idx = 0;

    ifstream matrix_file(filenames["bases"]);
    vector<cx_mat> output;
    smatch match, dummy;
    vector< complex<double> > numbers;
    
    string content((istreambuf_iterator<char>(matrix_file)), 
        istreambuf_iterator<char>());
    
    // Search the configured basis_file for basis matrix elements
    while(regex_search(content, match, elem_regex)){

        // for(auto i = 0u; i < match.size(); ++i)
        //     cout << i << ": " << match[i].str() << endl;

        string z = match[1].str();
        // cout << z << endl;
        //if (count(z.begin(), z.end(), '.') > 1){

        // do the split between "whole" complex numbers and partial ones
        if (regex_search(z, dummy, complex_split)) {
            string a = match[3].str(), b = match[4].str();
            clean_input(a); clean_input(b);
            // cout << stod(a) << "+" << stod(b) << "i" << endl;
            numbers.push_back(complex<double>(stod(a),stod(b)));
        } else {
            clean_input(z);
            // if (isdigit(z[z.length()-1]) || (z[z.length()-1] == '.')) {
            if (isdigit(z[z.length()-1]) || (z[z.length()-1] == '.')) {
                // cout << stod(z) << endl;
                numbers.push_back(complex<double>(stod(z), 0.0)); 
            } else {
                // cout << stod(z) << "i" << endl;
                numbers.push_back(complex<double>(0.0, stod(z)));
            }
        }
        content = match.suffix().str();
        idx++;
    }

    if (m*t*k != (int)numbers.size()){
        log_msg("Failed to read the " + to_string(m*t*k) + " matrix elements (read " + to_string(numbers.size()) + \
            ") configured in '" + filenames["bases"] + "'!", "Error");
        // log_msg();
        exit(1);
    }

    // store the read matrix elements in a complex matrix datatype from Armadillo library
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

map<int, int> read_error_requirements(const string &filepath){
    map<int, int> snr_to_error;
    vector<int> snrs, errors;
    int line_num = 1;
    ifstream error_file(filenames["error_file"]);
    string line;
    string value;

    if (!error_file.good()) {
        log_msg("No error requirements file '" + filepath + "' found!", "Error");
        exit(0);
    }

    while(getline(error_file, line)){
        istringstream iss(line);
        while(getline(iss, value, ',')){
            clean_input(value);
            if (line_num == 1){
                snrs.push_back(strtol(value.c_str(), NULL, 10));
            } else if (line_num == 2){
                errors.push_back(strtol(value.c_str(), NULL, 10));
            } else break;
        }
        line_num++;
    }

    if (snrs.empty()){
        log_msg("Error requirements file '" + filepath + "' is empty!", "Error");
        exit(0);
    }

    int min = params["snr_min"];
    int max = params["snr_max"];
    int step = params["snr_step"];
    unsigned required_size = (max-min)/step + 1;

    if (snrs.size() != required_size){
        log_msg("Error requirements file '" + filepath + "' has incorrect amount of values!", "Error");
        exit(0);
    }

    if (snrs.size() != errors.size()){
        log_msg("Error requirements file '" + filepath + "' is invalid!", "Error");
        exit(0);
    }

    for (auto i = 0u; i < snrs.size(); i++){
        snr_to_error[snrs[i]] = errors[i];
    }

    return snr_to_error;
}