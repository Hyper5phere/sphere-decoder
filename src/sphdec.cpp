/*
 * Compile with: g++ sphdec.cpp -o sphdec -O2 -larmadillo -llapack -lblas -std=c++14
 */

#include <iostream>
#include <armadillo>
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>

using namespace std;
using namespace arma;

void clean_input(string &input){
    // clean out comments
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    // clean out white spaces
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
}

void create_config(string filepath){
/* Creates a default configuration file in case one does not exist */
    ofstream defconf(filepath);
    defconf << "basis_file=bases.txt        // Text file containing the basis matrices" << endl \
            << "word_file=words.txt         // Text file containing the words" << endl \
            << "output_file=output.txt      // Text file used for simulation output" << endl \
            << "code_size=4                 // Size of the code" << endl \
            << "code_length=8               // Length of the code" << endl \
            << "no_of_matrices=4            // Number of matrices" << endl \
            << "no_of_broadcast_antennas=5  // Number of broadcast antennas" << endl \
            << "no_of_receiver_antennas=5   // Number of receiver antennas"  << endl \
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
        create_config(filepath);
        config_file = ifstream(filepath);
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

int main(int argc, char** argv)
{
    string inputfile = "options.ini";
    if (argc > 1){
        inputfile = argv[1];
    }
    auto output = get_options(inputfile);
    for (auto const &item : output){
        cout << item.first << " = " << item.second << endl;
    }
    return 0;
}