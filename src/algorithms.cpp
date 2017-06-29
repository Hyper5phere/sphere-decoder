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
int sesd_sign(double x){
    return (x <= 0) ? -1 : 1;
}

/* Takes the squared Frobenius norm from a complex matrix A */
double frob_norm_squared(cx_mat A){
    double sum = 0;
    for (auto i = 0u; i < A.n_rows; i++)
        for (auto j = 0u; j < A.n_cols; j++)
            sum += norm(A(i,j)); 
    return sum;
}

/* Makes sure the upper triangular matrix R has only positive diagonal elements */
void process_qr(mat &Q, mat &R){
    for (auto i = 0u; i < R.n_cols; i++){
        if (R(i,i) < 0.0){
            // R.col(i) *= -1;
            // Q.row(i) *= -1;
            R.row(i) *= -1;
            Q.col(i) *= -1;
        }
    }
}

/* Vectorizes a complex matrix A to real vector (each row is concatenated) */
vec to_real_vector(cx_mat A){
    vec a(2*A.n_elem); 
    // as the imaginary parts are now independent real elements in the vector,
    // we need to double the number of elements in the returned vector
    int index = 0;
    for (auto i = 0u; i < A.n_rows; i++){
        for (auto j = 0u; j < A.n_cols; j++){
            auto z = A(i,j);
            a[index++] = z.real();
            a[index++] = z.imag();
        }
    }
    return a;
}

/* generates a random n x m complex matrix from uniform distribution */
cx_mat create_random_matrix(int n, int m, double mean, double variance){
    normal_distribution<double> distr(mean, sqrt(variance));
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
vector<pair<vector<int>,cx_mat>> create_codebook(const vector<cx_mat> &bases, int* symbolset){
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

    vector<pair<vector<int>,cx_mat>> codebook;
    // vector<cx_mat> codebook;
    vector<int> tmp;
    cx_mat X(m, t);
    
    if(samples > 0){
        int random_index = 0;
        uniform_int_distribution<int> dist(0,q-1);
        log_msg("Random sampled data vector combinations:");
        for (int j = 0; j < samples; j++){
            for (int i = 0; i < k; i++) {
                random_index = dist(mersenne_twister);
                tmp.push_back(symbolset[random_index]);
                X = X + symbolset[random_index]*bases[i];
            }
            log_msg(vec2str(tmp, tmp.size()));
            tmp.clear();
            codebook.push_back(make_pair(tmp, X));
            X.zeros();
        }
        cout << endl;
    } else {

        /* all possible combinations of code words */
        auto c = comb_wrapper(symbolset, k);

        log_msg("All possible data vector combinations:");
        for (const auto &symbols : c){
            for (int i = 0; i < k; i++)
                X = X + symbols[i]*bases[i];
            
            codebook.push_back(make_pair(symbols, X));
            X.zeros();

            log_msg(vec2str(symbols, k));
        }
        cout << endl;
    }
    return codebook;
}

/* Computes the average and maximum energy of given codebook X */
pair<double,double> code_energy(const vector<pair<vector<int>,cx_mat>> &X){
    double sum = 0, max = 0, tmp = 0, average = 0;
    int cs = (int) X.size();

    for (int i = 0; i < cs; i++){
        //tmp = pow(norm(X[i], "fro"), 2);
        tmp = frob_norm_squared(X[i].second);
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