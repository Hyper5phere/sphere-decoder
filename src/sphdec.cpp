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

/* local anonymous namespace for helper functions */
namespace {
    /* Decision feedback equalization on xt[i] */
    inline void step2(int i, int q, vec &xt, const vec &y, vec &delta, const vec &ksi, const mat &R) {
        xt[i] = round((y[i]-ksi[i])/R(i,i));
        if (xt[i] < 0) {
            xt[i] = 0;
            delta[i] = 1;
        } else if (xt[i] > q - 1) {
            xt[i] = q - 1;
            delta[i] = -1;
        } else {
            delta[i] = sesd_sign(y[i]-ksi[i]-R(i,i)*xt[i]);
        }
    }

	/* Schnorr-Euchner enumeration step for the sphere decoder */
	inline void step6(int i, vec &xt, vec &delta){
		xt[i] = xt[i] + delta[i];
		delta[i] = -delta[i] - sesd_sign(delta[i]);
	}

    /* Perform some basic checks to ensure that the sphere decoder works as intended */
    inline bool check(double radius, const mat &R, const vec &y){
        if (radius <= 0){
            log_msg("sphdec: negative squared initial radius given!", "Error");
            log_msg("aborting simulation round...");
            return false;
        }

        if (R.n_rows != R.n_cols){
            log_msg("sphdec: R is not a square matrix!", "Error");
            log_msg("aborting simulation round...");
            return false;
        }

        if (y.n_elem != R.n_cols){
            log_msg("sphdec: vector y dimension mismatch!", "Error");
            log_msg("aborting simulation round...");
            return false;
        }
        return true;
    }
}

/* Basic sphere decoder algorithm */
vector<int> sphdec(const vec &y, const mat &R, int &counter, double radius){ //, vector<cx_mat> bases){

    // Step 1
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1;

    vector<int> x(k); // found point
	vec xt(k), ksi(k), delta(k), dist(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros();
    counter = 0;
	
	double xidist = 0.0;
    bool found = false;

    if (!check(radius, R, y))
        return x;

    step2(i, q, xt, y, delta, ksi, R);

    while (true) {
        counter++;
        // Step 3.
        xidist = pow(y[i]-ksi[i]-R(i,i)*xt[i], 2);
        
        // cout << i << ": " << vec2str(xt, k) << ", xidist = " << xidist + dist[i] << ", C = " << radius << endl;

        if (radius < dist[i] + xidist) { // current point xt is outside the sphere
        	// Step 4.
        	if (i == k-1) {
        		break;
        	} else {
        		// Step 6.
                i++;
        		step6(i, xt, delta);
        	}
        } else { // we are inside the sphere
        	if (xt[i] < 0 || xt[i] > q - 1){ // we are outside the signal set boundaries
                if ((xt[i] < 0 && (xt[i] + delta[i]) > q - 1) || (xt[i] > q - 1 && (xt[i] + delta[i]) < 0)){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        step6(i, xt, delta);
                    }
                } else 
        	        step6(i, xt, delta);
            } else {
        		if (i > 0) {
                    ksi[i-1] = 0;
        			for (int j = i; j < k; j++)
        				ksi[i-1] += R(i-1, j)*xt[j];
        			dist[i-1] = dist[i] + xidist;
        			i--;
                    step2(i, q, xt, y, delta, ksi, R);
        		} else { // lattice point is found (Step 5)
        			radius = dist[0] + xidist;
                    found = true;
        			x = conv_to<vector<int>>::from(xt);
                    i++;
        			step6(i, xt, delta);
        		}
        	}
        }
    }
    if (found)
	    return x;
    else {
        log_msg("sphdec: point not found!", "Warning");
        return vector<int>(0);
    }
}

/* Sphere decoder algorithm that considers spherical shaping of the codebook */
vector<int> sphdec_spherical_shaping(const vec &y, const mat &HR, const mat &R, int &counter, double P, double radius){ //, vector<cx_mat> bases){

    // Step 1
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1;

    vector<int> x(k); // found point
    vec xt(k), ksi(k), delta(k), dist(k), curr(k), ener(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros(); curr.zeros(); ener.zeros();
    counter = 0;
    
    double xidist = 0.0, xiener = 0.0;
    bool found = false;

    if (!check(radius, HR, y))
        return x;

    step2(i, q, xt, y, delta, ksi, HR);

    while (true) {
        counter++;
        // Step 3.
        xidist = pow(y[i]-ksi[i]-HR(i,i)*xt[i], 2);
        
        // cout << i << ": " << vec2str(xt, k) << ", xidist = " << xidist + dist[i] << ", C = " << radius << endl;

        if (radius < dist[i] + xidist) { // current point xt is outside the sphere
            // Step 4.
            if (i == k-1) {
                break;
            } else {
                // Step 6.
                i++;
                step6(i, xt, delta);
            }
        } else { // we are inside the sphere

            xiener = pow((2*xt[i] - q + 1)*R(i,i) + curr[i], 2);

            // we are outside the signal set boundaries
            if (xt[i] < 0 || xt[i] > q - 1 || ener[i] + xiener > P + 10e-6) {
                if ((xt[i] < 0 && (xt[i] + delta[i]) > q - 1) || (xt[i] > q - 1 && (xt[i] + delta[i]) < 0)){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        step6(i, xt, delta);
                    }
                } else 
                    step6(i, xt, delta);
            } else {
                if (i > 0) {
                    ksi[i-1] = 0;
                    curr[i-1] = 0;
                    for (int j = i; j < k; j++){
                        ksi[i-1] += HR(i-1, j)*xt[j];
                        curr[i-1] += (2*xt[j] - q + 1)*R(i-1,j);
                    }
                    dist[i-1] = dist[i] + xidist;
                    ener[i-1] = ener[i] + xiener;

                    i--;
                    step2(i, q, xt, y, delta, ksi, HR);
                } else { // lattice point is found (Step 5)
                    radius = dist[0] + xidist;
                    found = true;
                    x = conv_to<vector<int>>::from(xt);
                    i++;
                    step6(i, xt, delta);
                }
            }
        }
    }
    if (found)
        return x;
    else {
        log_msg("sphdec: point not found!", "Warning");
        return vector<int>(0);
    }
}

/* Wrapper function for the sphere decoder to handle complex to real matrix conversion, 
   QR-decomposition and other mappings */
vector<int> sphdec_wrapper(const vector<cx_mat> &bases, const cx_mat basis_sum, const cx_mat &H, 
                           const cx_mat &X, const cx_mat &N, int &visited_nodes, double radius){

    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int P = params["spherical_shaping_max_power"];

    vector<int> x(k);
    cx_mat Y, Ynorm;                           // Helper complex matrices
    mat B(2*t*n, k), G(2*t*n,k), Q, R, Rorig;  // real matrices
    vec y, y2;                                 // helper and sphdec input vector
    
    Y = H*X + N;                            // Simulated code block that we would receive from MIMO-channel
    Ynorm = (Y + H*basis_sum*(q - 1))*0.5;  // normalize received matrix for the sphere decoder
    y = to_real_vector(Ynorm);              // convert Y to real vector

    for(int i = 0; i < k; i++){
        G.col(i) = to_real_vector(bases[i]);   // G = (X1 X2 ... Xk) = generator matrix of original lattice
        B.col(i) = to_real_vector(H*bases[i]); // B = (HX1 HX2 ... HXk) = generator matrix of faded lattice
    }

    qr_econ(Q, Rorig, G);  // QR-decomposition of G (omits zero rows in Rorig)
    process_qr(Q, Rorig);  // Make sure Rorig has positive diagonal elements
    Q.zeros();
    qr_econ(Q, R, B);  // QR-decomposition of B (omits zero rows in R)
    process_qr(Q, R);  // Make sure R has positive diagonal elements

    y2 = Q.st()*y;     // Map y to same basis as R

    if (P < 0)
        x = sphdec(y2, R, visited_nodes, radius); // Call the actual sphere decoder algorithm
    else
        x = sphdec_spherical_shaping(y2, R, Rorig, visited_nodes, P, radius);

    if (x.size() == 0) // point not found
        return vector<int>(0);

    for (int j = 0; j < k; j++)
        x[j] = 2*x[j] - q + 1;

    return x;
}