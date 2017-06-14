/*
 ===================================================================================================
 Name        : sphdec.cpp
 Author      : Pasi Pyrrö
 Version     : 1.0
 Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis
 Date        : 14.6.2017
 Description : Sphere Decoder in C++14
 Compilation : g++ sphdec.cpp -o sphdec -O2 -larmadillo -llapack -lblas -std=c++14
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
#include <algorithm>
#include <regex>
#include <locale>

#define NUM_OPTIONS 11 

using namespace std;
using namespace arma;

/* default filenames */
string options_filename = "settings.ini";
string basis_filename = "bases.txt";
string output_filename = "output.txt";

/* simulation parameters */
map<string, int> params;


void clean_input(string &input){
    // clean out comments
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    // clean out white spaces
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
}

void create_config(const string filepath){
/* Creates a default configuration file in case one does not exist */
    ofstream defconf(filepath);
    defconf << "// configuration settings and simulation parameters for the sphere decoder program //" << endl << endl \
            << "basis_file=bases.txt        // Text file containing the basis matrices" << endl \
            << "output_file=output.txt      // Text file used for simulation output" << endl \
            << "code_size=4                 // Size of the code" << endl \
            << "code_length=2               // Length of the code" << endl \
            << "no_of_matrices=2            // Number of basis matrices" << endl \
            << "no_of_broadcast_antennas=2  // Number of broadcast antennas" << endl \
            << "no_of_receiver_antennas=2   // Number of receiver antennas"  << endl \
            << "snr_min=6                   // Minimum value for signal-to-noise ratio" << endl \
            << "snr_max=12                  // Maximum value for signal-to-noise ratio" << endl \
            << "symbolset=3                 // List of admissable values for the signal vectors" << endl \
            << "required_errors=500         // Demand at minimum this many errors before the end of the simulation" << endl;
    defconf.close();
}

void configure(const string filepath) {
/* Read simulation options from a configuration file */
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
                    basis_filename = value;
                else if (key.compare("output_file") == 0)
                    output_filename = value;
                else {
                    if ((params[key] = strtol(value.c_str(), NULL, 10)) <= 0) {
                        cout << "[Error] Value for option '" << key << "' must be an positive integer!" << endl;
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
    if (lines < NUM_OPTIONS) {
        cout << "[Error] too few options spesified in the '" << filepath << "!" << endl;
        cout << "[Info] Consider deleting the default options file which will reset the program settings." << endl;
        exit(1);
    }

    return;
}

vector<cx_mat> read_matrices(){
/* reads k (m x n) matrices from the spesified basis_file */

    /* - regular expression pattern used to search 
     *   matrix elements in Wolfram Mathematica format
     * - ignores most white spaces and supports asterisk in front of 'I'
     * - pattern is also case insensitive
     * - requires the use of decimal marks even with integer coefficients
     */
    regex elem_regex(
        "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*\\*?[Ii]?)|"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*" // white space
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\*?[Ii]{1})\\s*([,}]{1})" // element ends in 'I' and comma or bracket
    );

    /* Some supported formats for the matrix elements:
     * 1) "re ,"
     * 2) "re + im *I,"
     * 3) "re - im *I,"
     * 4) "im *I," 
     */

    size_t m = params["no_of_receiver_antennas"];
    size_t n = params["code_length"];
    size_t k = params["no_of_matrices"];
    size_t idx = 0;

    ifstream matrix_file(basis_filename);
    vector<cx_mat> output;
    smatch match;
    vector< complex<double> > numbers;
    
    string content((istreambuf_iterator<char>(matrix_file)), 
        istreambuf_iterator<char>());
    
    // Search the configured basis_file for basis matrix elements
    while(regex_search(content, match, elem_regex)){

        // for(auto i = 0u; i < match.size(); ++i)
        //     cout << i << ": " << match[i].str() << endl;

        string z = match[1].str();
        if (count(z.begin(), z.end(), '.') > 1){
            string a = match[3].str(), b = match[4].str();
            clean_input(a); clean_input(b);
            // cout << stod(a) << "+" << stod(b) << "i" << endl;
            numbers.push_back(complex<double>(stod(a),stod(b)));
        } else {
            clean_input(z);
            if (isdigit(z[z.length()-1])) {
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

    if (m*n*k != numbers.size()){
        cout << "[Error] Failed to read configured " <<  m*n*k << " matrix elements!" << endl;
        exit(1);
    }

    // store the read matrix elements in a complex matrix datatype from Armadillo library
    idx = 0;
    cx_mat X(m,n);
    for (auto i = 0u; i < k; i++){    
        for (auto j = 0u; j < m; j++){
            for (auto s = 0u; s < n; s++){
                X(j,s) = numbers[idx];
                idx++;
            }
        }
        output.push_back(X);
        X = cx_mat(m, n);
    }
    
    return output;
}

pair<double,double> code_energy(const vector<cx_mat> &X){
    double sum = 0, max = 0, tmp = 0, average = 0;
    for (int i = 0; i < params["code_length"]; i++){
        tmp += norm(X[i], "fro");
        sum += tmp;
        if (tmp > max)
            max = tmp;
    }
    average = sum / params["code_size"];
    return make_pair(average, max);
}

int main(int argc, char** argv)
{
    string inputfile = options_filename;
    if (argc == 2){
        inputfile = argv[1];
    } else if (argc > 2) {
        cout << "Usage: " << argv[0] << " [settings_file*]" << endl;
        exit(0);
    }

    configure(inputfile);

    auto bases = read_matrices();
    for (auto const &base : bases){
        cout << base << endl;
    }

    auto e = code_energy(bases);
    cout << "Code Energy" << \
    endl << "-----------" << \
    endl << "Average: " << e.first << \
    endl << "Max: " << e.second << endl;

    return 0;
}
