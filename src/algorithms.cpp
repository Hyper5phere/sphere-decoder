/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : algorithms.cpp                                                                        *
 * Project     : Schnorr-Euchnerr sphere decoder simulation for space-time lattice codes               *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++11                                                                                 *
 * Description : All custom mathematical algorithms not found in Armadillo or STL are collected here   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ARMA_NO_DEBUG /* disable Armadillo bound checks for addiotional speed */
#define FPLLL_WITH_ZDOUBLE
#define FPLLL_WITH_LONG_DOUBLE

#include <iostream>
#include <armadillo> /* linear algebra library */
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include <set>
#include <random>
#include <omp.h>
#include <tuple>

/* library for LLL reduction (optional) */
#ifdef USE_LLL
#include <fplll.h>
#endif
/* required packages to install:
 *	 libmpfr-dev
 *   libgmp-dev
 */

#include "algorithms.hpp"
#include "misc.hpp"

using namespace std;
using namespace arma;

/* Signum function for schorr-euchnerr sphere decoder (sesd) 
 * ---------------------------------------------------------
 * Note that sesd_sign(0) = -1
 */
int sesd_sign(double x) {
    return (x <= 0) ? -1 : 1;
}


/* compute greatest common divisor of a and b
 * credits go to: https://codereview.stackexchange.com/a/66735
 */
int gcd(int a, int b) {
    return b == 0 ? a : gcd(b, a % b);
}

/* Rounds x to nearest integer in S
 * --------------------------------
 * E.g. let S = 4-PAM and x = 1.8 then 
 * nearest_symbol(1.8) = 3
 * cf. round(1.8) = 2
 */
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

/* A simple heuristic function based on the idea that the volume of the hypersphere divided by the volume of the fundamental region of the lattice 
 * should roughly equal to the number of lattice points inside the hypersphere
 * i.e. this function estimates the squared radius required for the hypersphere in order to have 2^s codewords (lattice points) inside it
 * (for q-PAM we need to multiply the generator matrix by 2 to get the correct volume for the fundamental region)
 */
double estimate_squared_radius(const mat &G, int s){
    int n = params["no_of_matrices"];
    double pi = 3.1415926535897;
    // cout << "lattice constant: " << to_string(sqrt(det(4*G.t()*G))) << endl;
    // cout << "volume of the sphere: " << to_string(pow(pi, n/2.0)*pow(dparams["spherical_shaping_max_power"], n/2)/tgamma(n/2.0 + 1.0)) << endl;

    /* The basic idea: 2^s = vol(Sphere(R))/vol(Lambda) 
     * --> Solve for R when s and Lambda (G) are given
     */

    /* This algorithm uses this equation for the volume of the n-ball: 
     * https://en.wikipedia.org/wiki/Volume_of_an_n-ball
     */
    return pow(pow(2.0, s)*tgamma(n/2.0 + 1.0)*sqrt(det(4*G.t()*G))*pow(pi, n/-2.0), 2.0/n);
}

/* Generates a new symbolset where elements are between given lower and upper bounds */
vector<int> slice_symbset(const vector<int> &S, double lb, double ub){
    vector<int> retval;
    for (const int q : S){
        if (lb <= (double)q && (double)q <= ub) {
            retval.push_back(q);
        }
    }
    return retval;
}

/* Takes the squared Frobenius norm from a complex matrix A */
double frob_norm_squared(const cx_mat &A){
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

/* Create the lattice generator matrix G out of basis matrices B_i
 * i.e. G = vec(B_i), i = 1,...,k
 */
cx_mat create_generator_matrix(const vector<cx_mat> &bases){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    cx_mat G(k, m*t);
    for(int i = 0; i < k; i++){
        G.row(i) = vectorise(bases[i], 1);
    }
    return G.st();
}


/***** DEPRECATED *****
 * Use 
 * G_real = to_real_matrix(create_generator_matrix(bases));
 * instead 
 **********************/

/* Create the real valued lattice generator matrix G out of basis matrices */
// mat create_real_generator_matrix(const vector<cx_mat> &bases){
//     int m = params["no_of_transmit_antennas"];
//     int t = params["time_slots"];
//     int k = params["no_of_matrices"];
//     mat G(2*m*t,k);
//     for(int i = 0; i < k; i++){    
//         G.col(i) = to_real_vector(bases[i]);
//     }
//     return G;
// }


/* Inverse of the create_generator_matrix() function 
 * Useful for generating basis matrix files in conjunction with
 * output_complex_matrix()
 * See main.cpp for an example of this use case
 */
vector<cx_mat> generator_to_bases(const cx_mat &G){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int idx = 0;

    vector<cx_mat> bases;
    cx_mat X(m, t);
    cx_vec col(m*t);

    for (auto i = 0u; i < G.n_cols; i++){
        X.zeros();
        col = G.col(i);
        for (int j = 0; j < m; j++){
            for (int s = 0; s < t; s++){
                X(j,s) = col[idx];
                idx++;
            }
        }
        idx = 0;
        bases.push_back(X);
    }
    return bases;
}

/* Vectorizes a complex matrix A to real vector (each row is concatenated) */
vec to_real_vector(const cx_mat &A, bool row_wise){
    vec a(2*A.n_elem); 
    // as the imaginary parts are now independent real elements in the vector,
    // we need to double the number of elements in the returned vector
    int index = 0;
    if (row_wise){
        for (auto i = 0u; i < A.n_rows; i++){
            for (auto j = 0u; j < A.n_cols; j++){
                auto z = A(i,j);
                a[index++] = z.real();
                a[index++] = z.imag();
            }
        }
    } else { /* concatenate columns instead */
        for (auto i = 0u; i < A.n_cols; i++){
            for (auto j = 0u; j < A.n_rows; j++){
                auto z = A(j,i);
                a[index++] = z.real();
                a[index++] = z.imag();
            }
        }
    }
    return a;
}

/* Converts a complex matrix A to its real representation */
mat to_real_matrix(const cx_mat &A){
    mat B(2*A.n_rows, A.n_cols);
    for(auto i = 0u; i < A.n_cols; i++){   
        B.col(i) = to_real_vector(A.col(i));
    }
    return B;
}

/* Converts a real matrix constructed out of a complex matrix back to complex */
cx_mat to_complex_matrix(const mat &A){
    if (A.n_rows % 2 != 0){
        log_msg("to_complex_matrix: the input matrix needs to have an even number of rows!", "Error");
        exit(1);
    }

    cx_mat B(A.n_rows/2, A.n_cols);
    // int index = 0;
    for (auto j = 0u; j < A.n_cols; j++){
        for (auto i = 0u; i < B.n_rows; i++){
            B(i,j) = complex<double>(A(2*i,j), A(2*i+1,j));
        }
        // index = 0;
    }
    return B;
}

/* generates a random n x m 'full' complex matrix from normal distribution */
cx_mat create_random_matrix(int n, int m, double mean, double variance){
    normal_distribution<double> distr(mean, sqrt(variance));
    cx_mat A(n,m);
    return A.imbue([&]() { /* Armadillo magic with lambda function */
        return complex<double>(distr(mersenne_twister), distr(mersenne_twister));
    });
}

/* generates a random n x n diagonal complex matrix from normal distribution 
 * Used for generating the SISO channel matrix
 */
cx_mat create_random_diag_matrix(int n, double mean, double variance){
    normal_distribution<double> distr(mean, sqrt(variance));
    cx_vec diag(n);
    cx_mat A(n,n);

    A.zeros();
    diag.imbue([&]() { /* Armadillo magic with lambda function */
        return complex<double>(distr(mersenne_twister), 0.0);
    });
    A.diag() = diag;
    return A;
}

/* Creates q-PAM symbol set (integer vector) 
 * q is usually some power of two, i.e. q = 2, 4, 8, 16...
 * Example: 4-PAM = {-3, -1, 1, 3}
 */
vector<int> create_symbolset(int q){
    vector<int> symbset(q);
    for (int u = 0; u < q; u++) {
        symbset[u] = 2*u - q + 1;
    }
    return symbset; 
}

/* Calculates all possible combinations of elements for code vector _a_
   given the set of feasible symbols (element values)
   Used for generating the whole codebook, rarely needed though */
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

/* Attempts to take account the probability bias when picking vector coefficients uniform randomly 
   from hyperplane intervals within a hypersphere */
vector<int> create_unbiased_subset(const mat &R, const vector<int> &subset, vec xt,
                                   vec ener, vec curr, double radius, int i){
    if (i <= 0) return subset; /* can't perform lookahead for the lowest dimension */

    vector<int> unbiased_subset; //, delta_vec;
    unbiased_subset.reserve(4*subset.size());

    double xiener = 0.0;
    int ub = 0, lb = 0, delta = 0; //, delta_sum = 0, bias_correction_amount = 0;
    int k = params["no_of_matrices"];

    // cout << i << endl;
    // cout << vec2str(xt, xt.size()) << endl;
    // cout << vec2str(ener, ener.size()) << endl;
    // cout << vec2str(curr, curr.size()) << endl << endl;

    for (const int elem : subset) {
        xt[i] = elem;
        xiener = pow(R(i,i)*xt[i] + curr[i], 2);

        for (int j = i; j < k; j++)
            curr[i-1] += xt[j]*R(i-1,j);
        ener[i-1] = ener[i] + xiener;

        lb = nearest_symbol(-(sqrt(radius - ener[i-1]) + curr[i-1])/R(i-1,i-1), subset);
        ub = nearest_symbol((sqrt(radius - ener[i-1]) - curr[i-1])/R(i-1,i-1), subset);
        delta = (ub - lb)/2;

        // cout << lb << ", " << ub << ", " << delta << endl;

        for (int d = 0; d < delta; d++)
            unbiased_subset.push_back(elem);
        // delta_vec.push_back(delta);
        // delta_sum += delta;
    }
    // cout << vec2str(xt, xt.size()) << endl;
    // cout << vec2str(ener, ener.size()) << endl;
    // cout << vec2str(curr, curr.size()) << endl << endl;
    // cout << endl;

    // for (auto j = 0u; j < subset.size(); j++){
    //     bias_correction_amount = delta_sum/gcd(delta_sum, delta_vec[j]);
    //     for (int s = 0; s < bias_correction_amount; s++)
    //         unbiased_subset.push_back(subset[j]);
    // }
    // cout << vec2str(unbiased_subset, unbiased_subset.size()) << endl << endl;
    return unbiased_subset;
}

/* Creates a random codeword from basis matrices B_i and symbolset x-PAM within given radius (spherical shaping) */
pair<vector<int>, cx_mat> create_random_spherical_codeword(const vector<cx_mat> &bases, const mat &R, const vector<int> &S, double radius){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];

    cx_mat X(m, t);
    vector<int> coeffs(k), subset;
    X.zeros();

    vec xt(k), curr(k), ener(k);
    xt.zeros(); curr.zeros(); ener.zeros();
    int i = k-1;
    double xiener = 0.0;
    double lb = 0.0; // lower bound
    double ub = 0.0; // upper bound

    /* Randomly pick each component of the codeword and check if it's still within the energy bound in each dimension */
    while (i >= 0){
        lb = -(sqrt(radius - ener[i]) + curr[i])/R(i,i);
        ub = (sqrt(radius - ener[i]) - curr[i])/R(i,i);
        // uniform_real_distribution<double> xirange(lb, ub);
        // xt[i] = nearest_symbol(xirange(mersenne_twister), S); // probability bias fix this
        subset = slice_symbset(S, lb, ub);
        // if (i >= 3 /*max(k-13, 3)*/ && !subset.empty())
        //     subset = create_unbiased_subset(R, subset, xt, ener, curr, radius, i);

        if (subset.empty()) { /* No lattice points are within energy bounds in this dimension */
            i = k-1;
            xt.zeros(); curr.zeros(); ener.zeros();
            continue;
        }
        
        xt[i] = pick_uniform(subset);
        xiener = pow(R(i,i)*xt[i] + curr[i], 2);

        // if (xiener + ener[i] < radius + 1e-6) {  
        if (i > 0) {                     
            for (int j = i; j < k; j++)
                curr[i-1] += xt[j]*R(i-1,j);
            ener[i-1] = ener[i] + xiener;      
        }  
        i--;
        // } else { /* we're outside the energy bound, start over */
        //     log_msg("test!", "Error");
        //     i = k-1;
        //     xt.zeros(); curr.zeros(); ener.zeros();
        // }
    }
    coeffs = conv_to<vector<int>>::from(xt);
    for (int i = 0; i < k; i++) {
        X = X + coeffs[i]*bases[i];
    }
    return make_pair(coeffs, X);
}

/* Creates a codebook from basis matrices B_i and symbolset x-PAM within given radius (spherical shaping) */
// vector<pair<vector<int>, cx_mat>> create_spherical_codebook(const vector<cx_mat> &bases, const mat &R, const vector<int> &S, double radius){
vector<pair<vector<int>, cx_mat>> create_spherical_codebook(const vector<cx_mat> &bases, const mat &R, const vector<int> &S, 
															double radius, vec xt, int dim, double dist) {
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];

    cx_mat X(m, t, fill::zeros);
    vector<int> coeffs(k); //, subset;
    vector<pair<vector<int>,cx_mat>> codebook, part;

    // vec xt(k), curr(k), ener(k), indices(k);
    // xt.zeros(); curr.zeros(); ener.zeros(); indices.zeros();
    // int i = k-1;
    int i = dim-1;
    double xidist = 0;
    // double xiener = 0.0;
    // double lb = 0.0; // lower bound
    // double ub = 0.0; // upper bound
    // bool nextlevel = true;

    // while (true) {
    // 	if (nextlevel) {
	   //      lb = -(sqrt(radius - ener[i]) + curr[i])/R(i,i);
	   //      ub = (sqrt(radius - ener[i]) - curr[i])/R(i,i);
	   //      subset = slice_symbset(S, lb, ub);
	   //      if (subset.empty()) { /* No lattice points are within energy bounds in this dimension */
	   //      	cout << "empty interval!" << endl;
	   //      	if (i == k - 1) break;
	   //      	else {
	   //      		indices[i] = 0;
	   //      		// indices[i+1]++;
		  //           i++;
		  //           continue;
		  //       }
	   //      }
	   //      nextlevel = false;
	   //  }
    //   	xt[i] = subset[indices[i]];

    //     xiener = pow(R(i,i)*xt[i] + curr[i], 2);

    //     cout << i << ": " << vec2str(xt, k) << ", xiener = " << xiener + ener[i] << " interval: [" << lb << ", " << ub << "], xt[i] = " << xt[i] << endl;
    //     cout << "subset: " << vec2str(subset, subset.size()) << endl;

    //     // if (xiener + ener[i] < radius + 1e-6) {  
    //     if (i > 0) {
    //     	curr[i-1] = 0;                     
    //         for (int j = i; j < k; j++)
    //             curr[i-1] += xt[j]*R(i-1,j);
    //         ener[i-1] = ener[i] + xiener;  
    //         i--;
    //         nextlevel = true;
    //         cout << "lowering dimension!" << endl;
    //         // continue;    
    //     } else {
    //     	coeffs = conv_to<vector<int>>::from(xt);
    //     	X.zeros();
    //     	for (int i = 0; i < k; i++) {
		  //       X = X + coeffs[i]*bases[i];
		  //   }
		  //   cout << "found point!" << endl;
		  //   codebook.push_back(make_pair(coeffs, X));
    //     } 
    //     // }
    //     if (indices[i] >= subset.size() - 1){
    //     	cout << "subset clear!" << endl;
    //     	if (i == k - 1) break;
	   //  	else {
	   //  		indices[i] = 0;
	   //  		// indices[i+1]++;
	   //      	i++;
	   //      	nextlevel = true;
	   //      }
    //     } else indices[i]++;   	    
    // }
    // return codebook;

    // #pragma omp parallel for
    for (auto j = 0u; j < S.size(); j++) {
        vec tmp = xt;
        tmp[i] = S[j];
        xidist = pow(dot(R.row(i).subvec(i, k-1), tmp.subvec(i, k-1)), 2) + dist;
        if (xidist <= radius) {
            if (i > 0) {
            	part = create_spherical_codebook(bases, R, S, radius, tmp, dim-1, xidist);
                codebook.insert(codebook.end(), part.begin(), part.end());
            } else {
            	coeffs = conv_to<vector<int>>::from(tmp);
            	X.zeros();
	        	for (int s = 0; s < k; s++) {
			        X = X + coeffs[s]*bases[s];
			    }
			    // cout << vec2str(coeffs, k) << endl;
                codebook.push_back(make_pair(coeffs, X));
            }
        }
    }
    return codebook;
}


/* Creates a codebook (set of X matrices) from basis matrices B_i and symbolset q-PAM */
vector<pair<vector<int>,cx_mat>> create_codebook(const vector<cx_mat> &bases, const mat &R, const vector<int> &symbolset){
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int cs = params["codebook_size"];
    int samples = params["energy_estimation_samples"];
    double P = dparams["spherical_shaping_max_power"];

    vector<pair<vector<int>,cx_mat>> codebook;
    pair<vector<int>,cx_mat> code;
    vector<int> tmp;
    cx_mat X(m, t);
    X.zeros();
    
    if(samples > 0){
        for (int j = 0; j < samples; j++){
            if (P > 0)
                code = create_random_spherical_codeword(bases, R, symbolset, P);
            else
                code = create_random_codeword(bases, symbolset);
            codebook.push_back(code);
        }
    } else {
    	if (P > 0)
        	return create_spherical_codebook(bases, R, symbolset, P, vec(k, fill::zeros), k, 0);

        if (cs > 1e6)
            log_msg("create_codebook: Generating a codebook of size " + to_string(cs) + ", this can be really slow!", "Warning");

        /* all possible combinations of code words */
        auto c = comb_wrapper(symbolset, k);

        // log_msg("All possible data vector combinations:");
        for (const auto &symbols : c){
            X.zeros();
            for (int i = 0; i < k; i++)
                X = X + symbols[i]*bases[i];
            // log_msg(vec2str(symbols, symbols.size()));
            // if (frob_norm_squared(X) > P + 1e-6 && P > 0) continue;
            codebook.push_back(make_pair(symbols, X));
        }
    }
    return codebook;
}

/* Computes the average and maximum energy of given codebook X */
pair<double,double> code_energy(const vector<pair<vector<int>,cx_mat>> &X){
    double sum = 0, max = 0, tmp = 0, average = 0;
    int cs = (int) X.size();

    for (int i = 0; i < cs; i++){
        tmp = frob_norm_squared(X[i].second);
        sum += tmp;
        if (tmp > max)
            max = tmp;
    }
    average = sum / cs;
    return make_pair(average, max);
}

/* Returns the shortest lattice basis vector from the lattice generator matrix G*/
cx_vec shortest_basis_vector(const cx_mat &G){
    cx_vec shortest(G.n_rows, fill::zeros);
    double len = 10e9;
    for (auto i = 0u; i < G.n_cols; i++){
        if (len > norm(G.col(i), 2)) {
            shortest = G.col(i);
            len = norm(shortest, 2);
        }
    }
    return shortest;
}

/* Check if x belongs to the right coset (when using coset encoding)
 * i.e. check if the difference vector sent_x - decoded_x belongs in 
 * the given sublattice of the codebook carved out of lattice G_b
 */
bool coset_check(const cx_mat &Gb, const cx_mat &invGe, const Col<int> diff){
    cx_vec x = Gb*diff;
    cx_vec lambda = invGe*x;

    vec lambda_real = to_real_vector(lambda);
    
    /* check if lambda_real has only integer entries */
    for (const auto &l : lambda_real)
        if (fabs(l - round(l)) > 10e-5) /* tolerance */
            return false;

    return true;
}

/* calculates different rates related to the coset encoded block code */
tuple<double, double, double> code_rates(const cx_mat &Gb, const cx_mat &Ge){
    mat RGb = to_real_matrix(Gb);
    mat RGe = to_real_matrix(Ge);
    // cout << "Vol(Ge): " << sqrt(det(RGe.t()*RGe)) << endl;
    // cout << "Vol(Gb): " << sqrt(det(RGb.t()*RGb)) << endl;
    double r = log2(params["codebook_size"]);                          // overall rate
    double rt = log2(sqrt(det(RGe.t()*RGe))/sqrt(det(RGb.t()*RGb)));   // transmission rate
    double rc = r - rt;                                                // confusion rate
    return make_tuple(r, rt, rc); 
}

/* algorithm for counting the lattice points inside a _dim_ dimensional hypersphere of radius _radius_ 
 * considers only a single radius, but probably faster than the function below thanks to parallelisation
 */
int count_points(const mat &R, const vector<int> &S, double radius, vec xt, int dim, double dist){

    int k = params["no_of_matrices"];
    int i = dim - 1;

    vector<int> counters(omp_get_max_threads());

    #pragma omp parallel for
    for (auto j = 0u; j < S.size(); j++) {
        vec tmp = xt;
        tmp[i] = S[j];
        double xidist = pow(dot(R.row(i).subvec(i, k-1), tmp.subvec(i, k-1)), 2) + dist;
        if (xidist <= radius) {
            if (i > 0) {
                counters[omp_get_thread_num()] += count_points(R, S, radius, tmp, dim-1, xidist);
            } else {
                counters[omp_get_thread_num()]++;
            }
        }
    }
    return accumulate(counters.begin(), counters.end(), 0);
}


/* algorithm for counting the lattice points inside a _dim_ dimensional hypersphere of radiuses r_1, r_2, ..., r_n */
vector<int> count_points_many_radiuses(const mat &R, const vector<int> &S, vector<double> radiuses, vec xt, int dim, double dist){

    int k = params["no_of_matrices"];
    int i = dim - 1;

    vector<int> counters(radiuses.size());

    for (auto j = 0u; j < S.size(); j++){
        vec tmp = xt;
        tmp[i] = S[j];
        double xidist = pow(dot(R.row(i).subvec(i, k-1), tmp.subvec(i, k-1)), 2) + dist;
        
        if (xidist <= radiuses[0]){
            if (i > 0){
                vector<int> counts = count_points_many_radiuses(R, S, radiuses, tmp, dim-1, xidist);
                for (auto s = 0u; s < radiuses.size(); s++)
                    counters[s] += counts[s];
            } else {
                for (auto k = 0u; k < radiuses.size(); k++){
                    if (xidist <= radiuses[k])
                        counters[k]++;
                }
            }
        }  
    }
    return counters;
}

/* 
 * Should return the LLL reduced (as close to orthogonal as possible) basis for lattice generated by G 
 */
cx_mat LLL_reduction(const cx_mat &G) {
    #ifdef USE_LLL

    	mat G_real = to_real_matrix(G);

    	// ZZ_mat<double> fpG(G_real.n_rows, G_real.n_cols);
        // ZZ_mat<double> fpG(G_real.n_cols, G_real.n_rows);
        FP_mat<double> fpG(G_real.n_cols, G_real.n_rows);
    	for (auto i = 0u; i < G_real.n_rows; i++) {
    		for (auto j = 0u; j < G_real.n_cols; j++) {
    			// fpG(i, j) = G_real(i, j);
                fpG(j, i) = G_real(i, j);
    		}
    	}
        cout << fpG << endl << endl;
    	int status = lll_reduction(fpG, LLL_DEF_DELTA, LLL_DEF_ETA, LM_PROVED, FT_DEFAULT, 0, LLL_DEFAULT);
    	if (status != RED_SUCCESS) {
    	    log_msg("LLL reduction failed!", "Error"); //with error '" + to_string(get_red_status_str(status)) + "'", "Error");
    	    exit(1);
    	}
        cout << fpG << endl << endl;
    	for (auto i = 0u; i < G_real.n_rows; i++) {
    		for (auto j = 0u; j < G_real.n_cols; j++) {
    			// G_real(i, j) = fpG(i, j).get_d();
                G_real(i, j) = fpG(j, i).get_d();
    		}
    	}
        cout << G_real << endl;
    	return to_complex_matrix(G_real);

    #else /* Do nothing */

        return G;

    #endif
}