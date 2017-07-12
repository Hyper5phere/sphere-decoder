#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <string>
#include <vector>
#include <random>
#include <set>

#include "misc.hpp"

extern std::mt19937_64 mersenne_twister;

int sesd_sign(double x);
double frob_norm_squared(arma::cx_mat A);
void process_qr(arma::mat &Q, arma::mat &R);
arma::vec to_real_vector(arma::cx_mat A);
arma::cx_mat create_random_matrix(int n, int m, double mean, double variance);
// int* create_symbolset(int q);
std::vector<int> create_symbolset(int q);
void combinations(parallel_set< std::vector<int> > &comblist, const std::vector<int> &symbset, std::vector<int> comb, int dim);
std::set< std::vector<int> > comb_wrapper(const std::vector<int> &symbset, int vector_len);
// std::vector<arma::cx_mat> create_codebook(const std::vector<arma::cx_mat> &bases, int* symbolset);
std::pair<std::vector<int>, arma::cx_mat> create_random_codeword(const std::vector<arma::cx_mat> &bases, const std::vector<int> &symbolset);
std::vector<std::pair<std::vector<int>,arma::cx_mat>> create_codebook(const std::vector<arma::cx_mat> &bases, const std::vector<int> &symbolset);
// std::pair<double,double> code_energy(const std::vector<arma::cx_mat> X);
std::pair<double,double> code_energy(const std::vector<std::pair<std::vector<int>,arma::cx_mat>> &X);


#endif /* ALGORITHMS_HPP */