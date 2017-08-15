#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <string>
#include <vector>
#include <random>
#include <set>
#include <tuple>

#include "misc.hpp"

extern std::mt19937_64 mersenne_twister;

int sesd_sign(double x);
double nearest_symbol(double x, const std::vector<int> &S);
double estimate_squared_radius(const arma::mat &G, int s);
double frob_norm_squared(const arma::cx_mat &A);
// double euclidean_norm(const std::vector<int> &x);
void process_qr(arma::mat &Q, arma::mat &R);
// arma::mat pseudo_inverse(arma::mat A);
arma::cx_mat create_generator_matrix(const std::vector<arma::cx_mat> &bases);
// arma::mat create_real_generator_matrix(const std::vector<arma::cx_mat> &bases);
std::vector<arma::cx_mat> generator_to_bases(const arma::cx_mat &G);
arma::vec to_real_vector(const arma::cx_mat &A, bool row_wise = true);
arma::mat to_real_matrix(const arma::cx_mat &A);
arma::cx_mat to_complex_matrix(const arma::mat &A);
arma::cx_mat create_random_matrix(int n, int m, double mean, double variance);
arma::cx_mat create_random_diag_matrix(int n, double mean, double variance);
// int* create_symbolset(int q);
std::vector<int> create_symbolset(int q);
void combinations(parallel_set< std::vector<int> > &comblist, const std::vector<int> &symbset, std::vector<int> comb, int dim);
std::set< std::vector<int> > comb_wrapper(const std::vector<int> &symbset, int vector_len);
// std::vector<arma::cx_mat> create_codebook(const std::vector<arma::cx_mat> &bases, int* symbolset);
std::pair<std::vector<int>, arma::cx_mat> create_random_codeword(const std::vector<arma::cx_mat> &bases, const std::vector<int> &symbolset);
std::pair<std::vector<int>, arma::cx_mat> create_random_spherical_codeword(const std::vector<arma::cx_mat> &bases, 
                                                                           const arma::mat &R, const std::vector<int> &S, double radius);
std::vector<std::pair<std::vector<int>,arma::cx_mat>> create_codebook(const std::vector<arma::cx_mat> &bases, const arma::mat &R, 
                                                                      const std::vector<int> &symbolset);
// std::pair<double,double> code_energy(const std::vector<arma::cx_mat> X);
std::pair<double,double> code_energy(const std::vector<std::pair<std::vector<int>,arma::cx_mat>> &X);
arma::cx_vec shortest_basis_vector(const arma::cx_mat &G);
bool coset_check(const arma::cx_mat &Gb, const arma::cx_mat &invGe, const arma::Col<int> diff);
std::tuple<double, double, double> code_rates(const arma::cx_mat &Gb, const arma::cx_mat &Ge);
int count_points(const arma::mat &R, const std::vector<int> &S, double radius, arma::vec xt, int dim, double dist);
std::vector<int> count_points_many_radiuses(const arma::mat &R, const std::vector<int> &S, std::vector<double> radiuses, arma::vec xt, int dim, double dist);
arma::cx_mat LLL_reduction(arma::cx_mat G);


#endif /* ALGORITHMS_HPP */