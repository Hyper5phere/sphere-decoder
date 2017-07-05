#ifndef SPHDEC_HPP
#define SPHDEC_HPP

#include <armadillo>
#include <vector>
#include <complex>

std::vector<int> sphdec(double radius, const arma::vec &y, const arma::mat &R, int &counter); //, std::vector<arma::cx_mat> bases);
std::vector<int> sphdec_wrapper(const std::vector<arma::cx_mat> &bases, const arma::cx_mat basis_sum,\
	const arma::cx_mat &H, const arma::cx_mat &X, const arma::cx_mat &N, double radius, int &visited_nodes);

#endif /* SPHDEC_HPP */