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

using namespace std;
using namespace arma;

string default_filename = "options.ini";

void clean_input(string &input){
    // clean out comments
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    // clean out white spaces
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
    // input.erase(remove(input.begin(), input.end(), '+'), input.end());
}

void create_config(string filepath){
/* Creates a default configuration file in case one does not exist */
    ofstream defconf(filepath);
    defconf << "basis_file=bases.txt        // Text file containing the basis matrices" << endl \
            << "word_file=words.txt         // Text file containing the words" << endl \
            << "output_file=output.txt      // Text file used for simulation output" << endl \
            << "code_size=4                 // Size of the code" << endl \
            << "code_length=2               // Length of the code" << endl \
            << "no_of_matrices=1            // Number of basis matrices" << endl \
            << "no_of_broadcast_antennas=2  // Number of broadcast antennas" << endl \
            << "no_of_receiver_antennas=2   // Number of receiver antennas"  << endl \
            << "snr_min=6                   // Minimum value for signal-to-noise ratio" << endl \
            << "snr_max=12                  // Maximum value for signal-to-noise ratio" << endl \
            << "initial_radius=-1           // Initial search radius for the sphere decoder" << endl \
            << "symbolset=-3,-1,1,3         // List of admissable values for the signal vectors" << endl \
            << "required_errors=500         // Demand at minimum this many errors before the end of the simulation" << endl;
    defconf.close();
}

map<string, string> get_options(string filepath){
/* Read simulation options from a configuration file */
    map<string, string> output;
    ifstream config_file(filepath);

    if (!config_file.good()){
        cout << "Warning: no config file found, creating a default one..." << endl;
        create_config(default_filename);
        config_file = ifstream(default_filename);
    }

    string line;
    while (getline(config_file, line))
    {
        istringstream iss(line);
        string key;

        if(getline(iss, key, '=') )
        {
            string value;
            if(getline(iss, value)) {
                clean_input(value);
                if (!value.empty())
                    output[key] = value;
            } else {
                cout << "Warning: value for option '" << key << "' not spesified!" << endl;
            }
        }
    }
    return output;
}

vector<cx_mat> read_matrices(map<string, string> &options){

    regex elem_regex(
        "(([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*[Ii]?)|"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "\\s*"
        "([+-]?\\s*\\d*\\.?\\d+|[+-]?\\s*\\d+\\.?\\d*)" // parse any float
        "[Ii]{1})\\s*([,}]{1})"  
    );

    uint m = atoi(options.at("no_of_receiver_antennas").c_str());
    uint n = atoi(options.at("code_length").c_str());
    uint k = atoi(options.at("no_of_matrices").c_str());

    ifstream matrix_file(options.at("basis_file"));
    vector<cx_mat> output;
    smatch match;
    vector< complex<double> > numbers;
    
    string content((istreambuf_iterator<char>(matrix_file)), 
        istreambuf_iterator<char>());
    
    uint idx = 0;
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
        cout << "Failed to read configured " <<  m*n*k << " matrix elements!" << endl;
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
    string inputfile = default_filename;
    if (argc > 1){
        inputfile = argv[1];
    }
    auto output = get_options(inputfile);
    // for (auto const &item : output){
    //     cout << item.first << " = " << item.second << endl;
    // }

    auto bases = read_matrices(output);
    for (auto const &base : bases){
        cout << base << endl;
    }
    return 0;
}
