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
    // remove comments
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    // remove white spaces
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
}

/* Creates a default configuration file in case one does not exist */
void create_config(const string filepath){

    ofstream defconf(filepath);
    defconf << "// configuration settings and simulation parameters for the sphere decoder program //" << endl << endl \
            << "basis_file=bases.txt          // Text file containing the basis matrices" << endl \
            << "output_file=output.txt        // Text file used for simulation output" << endl \
            << "x-PAM=4                       // The size of the PAM signaling set" << endl \
            << "code_length=2                 // Dimension of the code words" << endl \
            << "energy_estimation_samples=-1  // Number of samples to make the code energy estimation (-1 = sample all)" << endl \
            << "no_of_matrices=2              // Number of basis matrices" << endl \
            << "no_of_broadcast_antennas=2    // Number of broadcast antennas" << endl \
            << "no_of_receiver_antennas=2     // Number of receiver antennas"  << endl \
            << "snr_min=6                     // Minimum value for signal-to-noise ratio" << endl \
            << "snr_max=12                    // Maximum value for signal-to-noise ratio" << endl \
            << "required_errors=500           // Demand at minimum this many errors before the simulation ends" << endl;
    defconf.close();
            // This is probably unnecessary
            // << "code_size=4                   // Size of the set containing all possible data vectors" << endl           
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
                    basis_filename = value;
                else if (key.compare("output_file") == 0)
                    output_filename = value;
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

    return;
}

/* reads k (m x n) matrices from the spesified basis_file */
vector<cx_mat> read_matrices(){

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

/* Creates a C-style integer array representation of an x-PAM symbolset */
int* create_symbolset(int q = params["x-PAM"]){
    int *symbset = (int*) malloc(q * sizeof(int));
    for (int u = 0; u < q; u++) {
        symbset[u] = 2*u - q + 1;
    }
    return symbset; // must be free()'d after use
}

/* recursively calculate all possible combinations of given symbolset
and return a list of lists of them */
// vector<vector<int> > combinations(int* symbolset, vector<int> comb, int dim){
    
//     if (dim == 0) {
//         vector<vector<int> > curr(1, vector<int>(params["x-PAM"]));
//         curr[0] = comb;
//         return curr;
//     }
//     else {
//         vector<vector<int> > curr(params["x-PAM"], vector<int>(params["x-PAM"]));
//         for (int i = 0; i < params["x-PAM"]; i++){
//             comb[dim] = symbolset[i];
//             curr[0] = comb;
//             auto tmp = combinations(symbolset, comb, dim-1);
//             for (int j = 0; j < params["x-PAM"]; j++)
//                 curr[j] = 
//             rest.insert(rest.end(), curr.begin(), curr.end());
//         }
//         return curr;
//     }
// }

/* Creates a codebook (set of X matrices) from basis matrices B_i and symbolset x-PAM */
vector<cx_mat> create_codebook(const vector<cx_mat> &bases, int* symbolset){
    // int q = params["x-PAM"];
    int m = params["no_of_receiver_antennas"];
    int n = params["code_length"];
    int k = params["no_of_matrices"];

    /* all possible combinations of code words */
    // TODO: Make this work somehow
    // auto combs = combinations(symbolset, vector<int>(symbolset[0]), n);

    /* lattice generator matrix G */
    cx_mat G(m*n,k);
    for(int i = 0; i < k; i++){
        G.col(i) = vectorise(bases[i]);
    }

    vector<cx_mat> codebook(k);
    // for (int j = 0; j < k; j++){
    //     cx_mat X(m,n);

    //     for (uint i = 0; i < combs.size(); i++){
    //         X = G*vec(combs[i]);
    //         // X = X + symbolset[i]*bases[i];
    //     }
    //     codebook[j] = X;
    // }
    return codebook;
}

/* Computes the average and maximum energy of given codebook X */
pair<double,double> code_energy(const vector<cx_mat> &X){
    double sum = 0, max = 0, tmp = 0, average = 0;

    for (uint i = 0; i < X.size(); i++){
        tmp += norm(X[i], "fro");
        sum += tmp;
        if (tmp > max)
            max = tmp;
    }

    average = sum / X.size();
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


    free(symbset);
    return 0;
}
