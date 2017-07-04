#ifndef SPHDEC_HPP
#define SPHDEC_HPP

#include <armadillo>
#include <vector>
#include <complex>

std::vector<int> sphdec(double radius, arma::vec &y, arma::mat &R, int &counter); //, std::vector<arma::cx_mat> bases);

#endif /* SPHDEC_HPP */