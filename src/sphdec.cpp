#include <iostream>
#include <armadillo> // linear algebra library
#include <vector>
#include <complex>
#include <string>

#include "sphdec.hpp"
#include "misc.hpp"
#include "algorithms.hpp"

using namespace std;
using namespace arma;

/* Sphere decoder algorithm */
vec sphdec(double radius, vec y, mat R, vector<cx_mat> bases){

	// int m = params["no_of_transmit_antennas"];
 //    int n = params["no_of_receiver_antennas"];
 //    int t = params["time_slots"];
 //    int k = params["no_of_matrices"];
    int q = params["x-PAM"];

	vec x(q); // coefficients for the lattice point we're looking for

	if (radius <= 0){
        log_msg("sphdec: negative squared initial radius given!", "Error");
        log_msg("aborting simulation round...");
        return x;
    }

    if (R.n_rows != R.n_cols){
    	log_msg("sphdec: R is not a square matrix!", "Error");
    	log_msg("aborting simulation round...");
    	return x;
    }

    if (y.n_elem != R.n_cols){
    	log_msg("sphdec: vector y dimension mismatch!", "Error");
    	log_msg("aborting simulation round...");
    	return x;
    }

	return x;
}