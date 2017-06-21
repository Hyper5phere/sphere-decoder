#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <string>
#include <set>
#include <random>

#include "algorithms.hpp"
#include "misc.hpp"

using namespace std;
using namespace arma;

/* signum function */
int sign(double x){
    return (x < 0) ? -1 : 1;
}

/* Takes the squared Frobenius norm from a complex matrix A */
double frob_norm_squared(cx_mat A){
    double sum = 0;
    for (auto i = 0u; i < A.n_rows; i++)
        for (auto j = 0u; j < A.n_cols; j++)
            sum += norm(A(i,j)); 
    return sum;
}

/* generates a random n x m complex matrix from uniform distribution */
cx_mat create_random_matrix(int n, int m, double mean, double variance){
    uniform_real_distribution<double> distr(mean, variance);
    cx_mat A(n,m);
    return A.imbue([&]() {
        return complex<double>(distr(mersenne_twister), distr(mersenne_twister));
    });
}

/* Creates a C-style integer array representation of an x-PAM symbolset */
int* create_symbolset(int q){
    int *symbset = (int*) malloc(q * sizeof(int));
    for (int u = 0; u < q; u++) {
        symbset[u] = 2*u - q + 1;
    }
    return symbset; // must be free()'d after use
}

/* Calculates all combinations of elements for code vector _a_
   given the set of feasible symbols (element values) */
void combinations(set< vector<int> > &comblist, vector<int> symbset, vector<int> comb, int dim){
    comblist.insert(comb);
    if (dim >= 0){
        for (const int symbol : symbset){
            comb[dim] = symbol;
            combinations(comblist, symbset, comb, dim-1);
        }
    }
}

/* Helper function for above combinations algorithm */
set< vector<int> > comb_wrapper(int* symbset, int vector_len){
    set< vector<int> > comblist;
    vector<int> init(vector_len);
    vector<int> symbset_v(symbset, symbset + params["x-PAM"]);
    for (int i = 0; i < vector_len; i++)
        init[i] = symbset[0];
    combinations(comblist, symbset_v, init, vector_len-1);
    return comblist;
}

/* Creates a codebook (set of X matrices) from basis matrices B_i and symbolset x-PAM */
vector<cx_mat> create_codebook(const vector<cx_mat> &bases, int* symbolset){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    // int cs = params["codebook_size"];
    int samples = params["energy_estimation_samples"];

    /* lattice generator matrix G (alternative approach) */
    // cx_mat G(m*n,k);
    // for(int i = 0; i < k; i++){
    //     G.col(i) = vectorise(bases[i]);
    // }

    vector<cx_mat> codebook;
    cx_mat X(m, t);
    
    if(samples > 0){
        int random_index = 0;
        uniform_int_distribution<int> dist(0,q-1);
        cout << "Random sampled data vector combinations:" << endl;
        for (int j = 0; j < samples; j++){
            cout << "{";
            for (int i = 0; i < k; i++) {
                random_index = dist(mersenne_twister);
                cout << symbolset[random_index] << " ";
                X = X + symbolset[random_index]*bases[i];
            }
            cout << "}" << endl;
            codebook.push_back(X);
            X.zeros();
            // cout << endl << X << endl;

            // cout << "{";
            // for (int s = 0; s < k - 1; s++)
            //     cout << symbolset[s] << ", ";
            // cout << symbolset[k-1] << "}" << endl;
        }
        cout << endl;
    } else {

        /* all possible combinations of code words */
        auto c = comb_wrapper(symbolset, k);

        cout << "All possible data vector combinations:" << endl;
        for (const auto &symbols : c){
            for (int i = 0; i < k; i++)
                X = X + symbols[i]*bases[i];
            
            codebook.push_back(X);
            X.zeros();

            cout << "{";
            for (int s = 0; s < k - 1; s++)
                cout << symbols[s] << ", ";
            cout << symbols[k-1] << "}" << endl;
        }
        cout << endl;
    }
    return codebook;
}

/* Computes the average and maximum energy of given codebook X */
pair<double,double> code_energy(const vector<cx_mat> X){
    double sum = 0, max = 0, tmp = 0, average = 0;
    int cs = (int) X.size();

    for (int i = 0; i < cs; i++){
        //tmp = pow(norm(X[i], "fro"), 2);
        tmp = frob_norm_squared(X[i]);
        sum += tmp;
        // cout << X[i] << endl;
        // cout << tmp << endl;
        // cout << sum << endl;
        if (tmp > max)
            max = tmp;
    }
    average = sum / cs;
    return make_pair(average, max);
}