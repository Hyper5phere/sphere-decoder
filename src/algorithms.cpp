#include <iostream>
#include <armadillo> // linear algebra library
#include <complex>
#include <vector>
#include <string>
#include <set>
#include <random>
#include <omp.h>


#include "algorithms.hpp"
#include "misc.hpp"

using namespace std;
using namespace arma;

/* signum function */
int sesd_sign(double x){
    return (x <= 0) ? -1 : 1;
}

/* Rounds x to nearest integer in S */
double nearest_symbol(double x, const vector<int> &S){
    double min = 10e6;
    double d = 0.0;
    double nearest = 0.0;

    for (const int symbol : S) {
        d = fabs(x - symbol);
        if (d < min) {
            nearest = symbol;
            min = d;
        }
    }
    return nearest;
}

int pick_uniform(const vector<int> S){
    uniform_int_distribution<int> dist(0, S.size()-1);
    return S[dist(mersenne_twister)];
}

vector<int> slice_symbset(const vector<int> &S, double lb, double ub){
    vector<int> retval;
    for (const int q : S){
        if (lb <= q && q <= ub) {
            retval.push_back(q);
        }
    }
    return retval;
}

/* Takes the squared Frobenius norm from a complex matrix A */
double frob_norm_squared(cx_mat A){
    double sum = 0;
    for (auto i = 0u; i < A.n_rows; i++)
        for (auto j = 0u; j < A.n_cols; j++)
            sum += norm(A(i,j)); 
    return sum;
}

double euclidean_norm(const vector<int> &x){
    double sum = 0;
    for (const int &elem : x)
        sum += pow(elem, 2);
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

// mat pseudo_inverse(mat A){
//     return inv(A.t()*A)*A.t();
// }

/* Create the lattice generator matrix G out of basis matrices */
cx_mat create_generator_matrix(const vector<cx_mat> &bases){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    cx_mat G(m*t,k);
    for(int i = 0; i < k; i++){
        G.col(i) = vectorise(bases[i]);  
    }
    return G;
}

/* Create the real valued lattice generator matrix G out of basis matrices */
mat create_real_generator_matrix(const vector<cx_mat> &bases){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    mat G(2*m*t,k);
    for(int i = 0; i < k; i++){    
        G.col(i) = to_real_vector(bases[i]);
    }
    return G;
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

/* generates a random n x n complex matrix from uniform distribution */
cx_mat create_random_diag_matrix(int n, double mean, double variance){
    normal_distribution<double> distr(mean, sqrt(variance));
    cx_vec diag(n);
    cx_mat A(n,n);

    A.zeros();
    diag.imbue([&]() {
        return complex<double>(distr(mersenne_twister), 0.0);
    });
    A.diag() = diag;
    return A;
}

// /* Creates a C-style integer array representation of an x-PAM symbolset */
// int* create_symbolset(int q){
//     int *symbset = (int*) malloc(q * sizeof(int));
//     for (int u = 0; u < q; u++) {
//         symbset[u] = 2*u - q + 1;
//     }
//     return symbset; // must be free()'d after use
// }

/* Creates q-PAM symbolset (integer array) */
vector<int> create_symbolset(int q){
    vector<int> symbset(q);
    for (int u = 0; u < q; u++) {
        symbset[u] = 2*u - q + 1;
    }
    return symbset; 
}

/* Calculates all combinations of elements for code vector _a_
   given the set of feasible symbols (element values) */
void combinations(parallel_set< vector<int> > &comblist, const vector<int> &symbset, vector<int> comb, int dim){
    comblist.par_insert(comb);
    if (dim >= 0){
        #pragma omp parallel 
        {
            vector<int> new_comb = comb;
            #pragma omp for
            for (size_t i = 0; i < symbset.size(); i++){
                new_comb[dim] = symbset[i];
                combinations(comblist, symbset, new_comb, dim-1);
            }
        }
    }
}

/* Helper function for above combinations algorithm */
set< vector<int> > comb_wrapper(const vector<int> &symbset, int vector_len){
    parallel_set< vector<int> > comblist;
    vector<int> init(vector_len);
    // vector<int> symbset_v(symbset, symbset + params["x-PAM"]);
    for (int i = 0; i < vector_len; i++)
        init[i] = symbset[0];
    combinations(comblist, symbset, init, vector_len-1);
    return comblist;
}

/* Creates a random codeword from basis matrices B_i and symbolset x-PAM */
pair<vector<int>, cx_mat> create_random_codeword(const vector<cx_mat> &bases, const vector<int> &symbolset){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int random_index = 0;

    cx_mat X(m, t);
    vector<int> coeffs(k);
    X.zeros();

    uniform_int_distribution<int> dist(0, q-1);

    for (int i = 0; i < k; i++) {
        random_index = dist(mersenne_twister);
        coeffs[i] = symbolset[random_index];
        X = X + symbolset[random_index]*bases[i];
    }

    return make_pair(coeffs, X);   
}

/* Creates a random codeword from basis matrices B_i and symbolset x-PAM within given radius (spherical shaping) */
pair<vector<int>, cx_mat> create_random_spherical_codeword(const vector<cx_mat> &bases, const mat &R, const vector<int> &S, double radius){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];

    cx_mat X(m, t);
    vector<int> coeffs(k);
    X.zeros();

    vec xt(k), curr(k), ener(k);
    xt.zeros(); curr.zeros(); ener.zeros();
    int i = k-1;
    double xiener = 0.0;
    double lb = 0.0; // lower bound
    double ub = 0.0; // upper bound

    while (i >= 0){
        lb = -(sqrt(radius - ener[i]) + curr[i])/R(i,i);
        ub = (sqrt(radius - ener[i]) - curr[i])/R(i,i);
        // uniform_real_distribution<double> xirange(lb, ub);
        // xt[i] = nearest_symbol(xirange(mersenne_twister), S); // probability fix this
        vector<int> subset = slice_symbset(S, lb, ub);
        if (subset.empty()){
            i = k-1;
            xt.zeros(); curr.zeros(); ener.zeros();
            continue;
        }
        xt[i] = pick_uniform(subset);
        xiener = pow(R(i,i)*xt[i] + curr[i], 2);

        if (xiener + ener[i] < radius + 1e-6) {  
            if (i > 0) {                     
                for (int j = i; j < k; j++)
                    curr[i-1] += xt[j]*R(i-1,j);
                ener[i-1] = ener[i] + xiener;      
            }  
            i--;
        } else {
            i = k-1;
            xt.zeros(); curr.zeros(); ener.zeros();
        }
    }
    coeffs = conv_to<vector<int>>::from(xt);
    for (int i = 0; i < k; i++) {
        X = X + coeffs[i]*bases[i];
    }
    // cout << vec2str(coeffs, coeffs.size()) << endl;
    return make_pair(coeffs, X);
}


/* Creates a codebook (set of X matrices) from basis matrices B_i and symbolset x-PAM */
vector<pair<vector<int>,cx_mat>> create_codebook(const vector<cx_mat> &bases, const mat &R, const vector<int> &symbolset){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    // int q = params["x-PAM"];
    int cs = params["codebook_size"];
    int samples = params["energy_estimation_samples"];
    double P = dparams["spherical_shaping_max_power"];

    /* lattice generator matrix G (alternative approach) */
    // cx_mat G(m*n,k);
    // for(int i = 0; i < k; i++){
    //     G.col(i) = vectorise(bases[i]);
    // }

    vector<pair<vector<int>,cx_mat>> codebook;
    pair<vector<int>,cx_mat> code;
    vector<int> tmp;
    cx_mat X(m, t);
    X.zeros();
    
    if(samples > 0){
        // int random_index = 0;
        // uniform_int_distribution<int> dist(0,q-1);
        // log_msg("Random sampled data vector combinations:");
        for (int j = 0; j < samples; j++){
            // for (int i = 0; i < k; i++) {
            //     random_index = dist(mersenne_twister);
            //     tmp.push_back(symbolset[random_index]);
            //     X = X + symbolset[random_index]*bases[i];
            // }
            // // log_msg(vec2str(tmp, tmp.size()));
            // tmp.clear();
            // do {
            //     code = create_random_codeword(bases, symbolset);
            // } while (/*euclidean_norm(code.first)*/ frob_norm_squared(code.second) > P + 1e-6 && P > 0); // do until we're inside spherical constellation
            if (P > 0)
                code = create_random_spherical_codeword(bases, R, symbolset, P);
            else
                code = create_random_codeword(bases, symbolset);
            codebook.push_back(code);
            // X.zeros();
        }
        // cout << endl;
    } else {
        if (cs > 1e6)
            log_msg("create_codebook: Generating a codebook of size " + to_string(cs) + "!", "Warning");

        /* all possible combinations of code words */
        auto c = comb_wrapper(symbolset, k);

        // log_msg("All possible data vector combinations:");
        for (const auto &symbols : c){
            // if (euclidean_norm(symbols) > P + 1e-6 && P > 0) continue; // We're outside the spherical constellation
            X.zeros();
            for (int i = 0; i < k; i++)
                X = X + symbols[i]*bases[i];

            // cout << "------------" << endl;
            // cout << X << endl;
            // log_msg(vec2str(symbols, symbols.size()));
            if (frob_norm_squared(X) > P + 1e-6 && P > 0) continue;
            codebook.push_back(make_pair(symbols, X));
            // cout << "added to codebook!" << endl;
    
        }
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

/* Check if x belongs to the right coset (when using coset encoding)
 * i.e. check if the difference vector sent_x - decoded_x belongs in 
 * the given sublattice of the codebook carved out of lattice G_b
 */
bool coset_check(const cx_mat &Gb, const cx_mat &invGe, const Col<int> diff){
    cx_vec x = Gb*diff;
    cx_vec lambda = invGe*x;

    vec lambda_real = to_real_vector(cx_mat(lambda));

    // cout << vec2str(lambda_real, lambda_real.size()) << endl;
    
    /* check if lambda_real has only integer entries */
    for (const auto &l : lambda_real)
        if (fabs(l - round(l)) > 10e-5) // tolerance
            return false;

    return true;
}


/* algorithm for counting the lattice points inside a _dim_ dimensional hypersphere of radius _radius_ */
int count_points(const mat &R, const vector<int> &S, double radius, vec xt, int dim, double dist){

    int k = params["no_of_matrices"];
    int i = dim - 1;

    vector<int> counters(omp_get_max_threads());

    #pragma omp parallel for
    for (auto j = 0u; j < S.size(); j++){
        vec tmp = xt;
        tmp[i] = S[j];
        double xidist = pow(dot(R.row(i).subvec(i, k-1), tmp.subvec(i, k-1)), 2) + dist;
        if (xidist <= radius){
            if (i > 0){
                counters[omp_get_thread_num()] += count_points(R, S, radius, tmp, dim-1, xidist);
            } else {
                // cout << vec2str(xt, xt.size()) << endl;
                counters[omp_get_thread_num()]++;
            }
        }
    }
    return accumulate(counters.begin(), counters.end(), 0);
}