// #define ARMA_NO_DEBUG // for speed

#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <regex>

#include "config.hpp"
#include "misc.hpp"

#define NUM_OPTIONS 12

using namespace std;
using namespace arma;

/* Creates a default configuration file in case one does not exist */
void create_config(const string filepath){
    // the file can contain comments similar to this line here
    ofstream defconf(filepath);
    defconf << "// configuration settings and simulation parameters for the sphere decoder program //" << endl << endl \
            << "basis_file=bases.txt            // Text file containing the basis matrices" << endl \
            << "output_file=output.txt          // Text file used for simulation output" << endl \
            << "x-PAM=4                         // The size of the PAM signaling set (even positive integer)" << endl \
            << "energy_estimation_samples=-1    // Number of samples to make the code energy estimation (-1 = sample all)" << endl \
            << "no_of_matrices=2                // Number of basis matrices (dimension of the data vectors)" << endl \
            << "time_slots=2                    // Number of time slots used in the code" << endl \
            << "no_of_transmit_antennas=2       // Number of transmit antennas" << endl \
            << "no_of_receiver_antennas=2       // Number of receiver antennas"  << endl \
            << "snr_min=6                       // Minimum value for signal-to-noise ratio" << endl \
            << "snr_max=12                      // Maximum value for signal-to-noise ratio" << endl \
            << "snr_step=2                      // Increase SNR by this value per each iteration" << endl \
            << "required_errors=500             // Demand at minimum this many errors before the simulation ends" << endl;
    defconf.close();     
}

/* Read simulation parameters from a settings file */
void configure(const string filepath) {
    map<string, int> output;
    ifstream config_file(filepath);

    if (!config_file.good() && filepath.compare(options_filename) != 0) {
        cout << "[Warning] No options file '" << filepath << "' found, using the default one..." << endl;
        config_file = ifstream(options_filename);
    } 

    if (!config_file.good()) {
        cout << "[Info] No default settings file found, creating a new one with default settings..." << endl;
        create_config(options_filename);
        cout << "[Info] Make your changes to the options file and rerun this program." << endl;
        cout << "[Info] Exiting..." << endl;
        exit(0);
    }
    
    string line;
    uint lines = 0;
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
            if (!value.empty()) {
                clean_input(value);
                if (key.compare("basis_file") == 0)
                    basis_filename = "../bases/" + value;
                else if (key.compare("output_file") == 0)
                    output_filename = "../output/" + value;
                else {
                    if ((params[key] = strtol(value.c_str(), NULL, 10)) == 0) {
                        cout << "[Error] Invalid value for option '" << key << "'" << endl;
                        exit(1);
                    }
                }
            } else {
                cout << "[Error] Value for option '" << key << "' not spesified!" << endl;
                exit(1);
            }
        }
        lines++;
    }
    if (lines < NUM_OPTIONS){
        cout << "[Error] Too few options spesified in the '" << filepath << "!" << endl;
        cout << "[Info] Consider deleting the default settings file which will reset the program settings." << endl;
        exit(1);
    } 
    if ((params["x-PAM"] <= 0) || (params["x-PAM"] % 2 == 1)){
        cout << "[Error] Invalid x-PAM signal set spesified!" << endl;
        cout << "[Info] x-PAM option accepts positive even integers." << endl;
        exit(1);
    }

    // helper variable (size of the code matrix set) calculated from the input parameters
    params["codebook_size"] = (int)pow(params["x-PAM"], params["no_of_matrices"]);
}

/* reads k (m x t) matrices from the spesified basis_file */
vector<cx_mat> read_matrices(){

    /* - regular expression pattern used to search 
     *   matrix elements in almost any known format (eg. Mathematica, Matlab)
     * - ignores most white spaces and supports asterisk in front of 'I'
     * - pattern is also case insensitive
     */
    regex elem_regex(
        "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*\\*?[Ii]?)|"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*" // white space
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\*?[Ii]{1})\\s*([,)};\\] \n]{1})" // element ends in comma, semicolon, newline, white space or closing bracket
    );

    /* Some supported formats for the matrix elements:
     * 1) "re ,"
     * 2) "re + im *I,"
     * 3) "re - im *I,"
     * 4) "im *I," 
     */

    // Used to determine whether the element has both real and complex part
    regex complex_split("\\d+\\.?\\s*[+-]{1}");

    size_t m = params["no_of_transmit_antennas"];
    size_t t = params["time_slots"];
    size_t k = params["no_of_matrices"];
    size_t idx = 0;

    ifstream matrix_file(basis_filename);
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

    if (m*t*k != numbers.size()){
        cout << "[Error] Failed to read the " <<  m*t*k << " matrix elements configured in '" << basis_filename << "'!" << endl;
        exit(1);
    }

    // store the read matrix elements in a complex matrix datatype from Armadillo library
    idx = 0;
    cx_mat X(m, t);
    for (auto i = 0u; i < k; i++){    
        for (auto j = 0u; j < m; j++){
            for (auto s = 0u; s < t; s++){
                X(j,s) = numbers[idx];
                idx++;
            }
        }
        output.push_back(X);
        X.zeros();
    }
    
    return output;
}