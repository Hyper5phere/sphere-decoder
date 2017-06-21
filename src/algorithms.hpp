#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <string>
#include <vector>
#include <random>
#include <set>

extern std::mt19937_64 mersenne_twister;

int sign(double x);
double frob_norm_squared(arma::cx_mat A);
arma::cx_mat create_random_matrix(int n, int m, double mean, double variance);
int* create_symbolset(int q);
void combinations(std::set< std::vector<int> > &comblist, std::vector<int> symbset, std::vector<int> comb, int dim);
std::set< std::vector<int> > comb_wrapper(int* symbset, int vector_len);
std::vector<arma::cx_mat> create_codebook(const std::vector<arma::cx_mat> &bases, int* symbolset);
std::pair<double,double> code_energy(const std::vector<arma::cx_mat> X);

#endif /* ALGORITHMS_HPP */