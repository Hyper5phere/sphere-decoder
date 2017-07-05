#ifndef SPHDEC_HPP
#define SPHDEC_HPP

#include <armadillo>
#include <vector>
#include <complex>

std::vector<int> sphdec(double radius, arma::vec &y, arma::mat &R, int &counter); //, std::vector<arma::cx_mat> bases);
std::vector<int> sphdec_wrapper(arma::cx_mat H, arma::cx_mat X, arma::cx_mat N, int &visited_nodes);

#endif /* SPHDEC_HPP */