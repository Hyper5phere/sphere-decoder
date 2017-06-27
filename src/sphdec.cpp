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

namespace {
	/* Schnorr-Euchner enumeration step for the sphere decoder */
	inline void step6(int i, vec &xt, vec &delta, bool &nextlevel){
		xt[i] = xt[i] + delta[i];
		delta[i] = -delta[i] - sign(delta[i]);
		nextlevel = false;
	}
}

/* Sphere decoder algorithm */
vec sphdec(double radius, vec y, mat R, vector<cx_mat> bases){

    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1;

	vec x(k), xt(k), ksi(k), delta(k), dist(k);
	
	double d = 0.0, ksitemp = 0.0;
	bool nextlevel = false;

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

    while (true) {
        cout << i << ": " << vec2str(xt, k) << endl;
    	// step 2.  
    	if (nextlevel) {
            xt[i] = (int)round((y[i]-ksi[i])/R(i,i));
            if (xt[i] < 0) {
                xt[i] = 0;
                delta[i] = 1;
            } else if (xt[i] > q - 1) {
                xt[i] = q - 1;
                delta[i] = -1;
            } else {
                delta[i] = sign((y[i]-ksi[i])-R(i,i)*xt[i]);
            }
        }

        // Step 3.
        d = pow(fabs(y[i]-ksi[i]-R(i,i)*xt[i]), 2);
        
        if (radius < dist[i] + d) { // current point xt is outside the sphere
        	// Step 4.
        	if (i == k-1) {
        		break;
        	} else {
        		// Step 6.
                i++;
        		step6(i, xt, delta, nextlevel);
        	}
        } else { // we are inside the sphere
        	if (xt[i] < 0 || xt[i] > q - 1){ // we are outside the signal set boundaries
                if ((xt[i] < 0 && (xt[i] + delta[i]) > q - 1) || (xt[i] > q - 1 && (xt[i] + delta[i]) < 0)){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        step6(i, xt, delta, nextlevel);
                    }
                } else 
        	        step6(i, xt, delta, nextlevel);
            } else {
        		if (i > 0) {
        			for (int j = i; j < k; j++)
        				ksitemp += R(i-1, j)*x[j];
        			ksi[i-1] = ksitemp;
        			dist[i-1] = dist[i] + d;
        			i--;
        			nextlevel = true;
        		} else { // lattice point is found (Step 5)
        			radius = dist[0] + d;
        			x = xt;
                    i++;
        			step6(i, xt, delta, nextlevel);
        		}
        	}
        }
        // break;
    }

	return x;
}