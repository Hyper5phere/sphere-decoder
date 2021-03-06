/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : sphdec.hpp                                                                            *
 * Project     : Planewalker - Schnorr-Euchner sphere decoder simulation for space-time lattice codes  *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++ (2011 or newer standard)                                                          *
 * Description : Function prototypes for sphdec.cpp                                                    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SPHDEC_HPP
#define SPHDEC_HPP

#include <armadillo>
#include <vector>
#include <complex>
#include <limits>

std::vector<int> sphdec(const arma::vec &y, const arma::mat &R, const std::vector<int> &S,
						int &counter, double radius=std::numeric_limits<double>::max());
std::vector<int> sphdec_spherical_shaping(const arma::vec &y, const arma::mat &HR, const arma::mat &R, const std::vector<int> &S,
										  int &counter, double P, double radius=std::numeric_limits<double>::max());
std::vector<int> sphdec_wrapper(const std::vector<arma::cx_mat> &bases, const arma::mat Rorig,
								const arma::cx_mat &H, const arma::cx_mat &X, const arma::cx_mat &N, 
								const std::vector<int> &symbset, int &visited_nodes, double radius=std::numeric_limits<double>::max());

#endif /* SPHDEC_HPP */