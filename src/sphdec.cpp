/*
 * Compile with: g++ sphdec.cpp -o sphdec -O2 -larmadillo -llapack -lblas -std=c++14
 */

// #define ARMA_NO_DEBUG // for speed

#include <iostream>
#include <armadillo>
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

string options_filename = "options.ini";
string basis_filename = "bases.txt";
string output_filename = "output.txt";


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
            << "no_of_matrices=1            // Number of basis matrices" << endl \
            << "no_of_broadcast_antennas=2  // Number of broadcast antennas" << endl \
            << "no_of_receiver_antennas=2   // Number of receiver antennas"  << endl \
            << "snr_min=6                   // Minimum value for signal-to-noise ratio" << endl \
            << "snr_max=12                  // Maximum value for signal-to-noise ratio" << endl \
            << "symbolset=3                 // List of admissable values for the signal vectors" << endl \
            << "required_errors=500         // Demand at minimum this many errors before the end of the simulation" << endl;
    defconf.close();
}

map<string, int> get_options(const string filepath){
/* Read simulation options from a configuration file */
    map<string, int> output;
    ifstream config_file(filepath);

    if (!config_file.good() && filepath.compare(options_filename) != 0){
        cout << "[Warning] No options file '" << filepath << "' found, using the default one..." << endl;
        config_file = ifstream(options_filename);
    } 

    if (!config_file.good()){
        cout << "[Info] No default options file found, creating a new one with default settings..." << endl;
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
                    if ((output[key] = strtol(value.c_str(), NULL, 10)) <= 0) {
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

    return output;
}

vector<cx_mat> read_matrices(const map<string, int> &options){
/* reads k nxm matrices from the spesified input text files */

    /* regular expression pattern used to search 
     * matrix elements in Wolfram Mathematica format
     */
    regex elem_regex(
        "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*[Ii]?)|"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*" // white space
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "[Ii]{1})\\s*([,}]{1})" // element ends in 'I' and comma or bracket
    );

    // uint m = atoi(options.at("no_of_receiver_antennas").c_str());
    // uint n = atoi(options.at("code_length").c_str());
    // uint k = atoi(options.at("no_of_matrices").c_str());
    size_t m = options.at("no_of_receiver_antennas");
    size_t n = options.at("code_length");
    size_t k = options.at("no_of_matrices");
    size_t idx = 0;

    // ifstream matrix_file(options.at("basis_file"));
    ifstream matrix_file(basis_filename);
    vector<cx_mat> output;
    smatch match;
    vector< complex<double> > numbers;
    
    string content((istreambuf_iterator<char>(matrix_file)), 
        istreambuf_iterator<char>());
    
    while(regex_search(content, match, elem_regex)){
        // cout << content << endl;
        if (!match.empty()) {
            // matches.push_back(match);
            // for(auto i = 0u; i < match.size(); ++i)
            //     cout << i << ": " << match[i].str() << endl;
            // cout << m[0].str() << endl;
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
        } else {
            return output;
        }
        content = match.suffix().str();
        idx++;
    }

    if (m*n*k != numbers.size()){
        cout << "[Error] Failed to read configured " <<  m*n*k << " matrix elements!" << endl;
        exit(1);
    }

    idx = 0;
    cx_mat X(m,n);
    for (auto i = 0u; i < k; i++){    
        // cout << i << X << endl;
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

int main(int argc, char** argv)
{
    string inputfile = options_filename;
    if (argc == 2){
        inputfile = argv[1];
    } else if (argc > 2) {
        cout << "Usage: " << argv[0] << " [options_file*]" << endl;
        exit(0);
    }
    auto opts = get_options(inputfile);
    // for (auto const &item : output){
    //     cout << item.first << " = " << item.second << endl;
    // }

    auto bases = read_matrices(opts);
    for (auto const &base : bases){
        cout << base << endl;
    }
    return 0;
}
