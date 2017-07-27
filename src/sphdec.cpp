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


/* Decision feedback equalization on xt[i] */
inline void step2(int i, int q, vec &xt, const vec &y, vec &delta, const vec &ksi, const mat &R, const vector<int> &S) {
    int spacing = abs(S[0] - S[1]);
    // xt[i] = round((y[i]-ksi[i])/R(i,i));
    xt[i] = nearest_symbol((y[i]-ksi[i])/R(i,i), S);
    if (xt[i] < S[0]) {
        xt[i] = S[0];
        delta[i] = spacing;
    } else if (xt[i] > S[q - 1]) {
        xt[i] = S[q - 1];
        delta[i] = -1*spacing;
    } else {
        delta[i] = spacing*sesd_sign(y[i]-ksi[i]-R(i,i)*xt[i]);
    }
}

/* Schnorr-Euchner enumeration step for the sphere decoder */
inline void step6(int i, vec &xt, vec &delta, const vector<int> &S){
    int spacing = abs(S[0] - S[1]);
	xt[i] = xt[i] + delta[i];
	delta[i] = -delta[i] - spacing*sesd_sign(delta[i]);
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


/* Basic sphere decoder algorithm */
vector<int> sphdec(const vec &y, const mat &R, const vector<int> &S, int &counter, double radius){

    // Step 1
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1;

    // cout << "sphdec called" << endl;

    vector<int> x(k); // candidate for the closest lattice point
	vec xt(k), ksi(k), delta(k), dist(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros();
    counter = 0;
	
	double xidist = 0.0;
    bool found = false;

    if (!check(radius, R, y))
        return x;

    step2(i, q, xt, y, delta, ksi, R, S);

    while (true) {

        if (exit_flag) break; // terminate simulations

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
        		step6(i, xt, delta, S);
        	}
        } else { // we are inside the sphere
        	if (xt[i] < S[0] || xt[i] > S[q - 1]){ // we are outside the signal set boundaries
                if ((xt[i] < S[0] && (xt[i] + delta[i]) > S[q - 1]) || (xt[i] > S[q - 1] && (xt[i] + delta[i]) < S[0])){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        step6(i, xt, delta, S);
                    }
                } else {
                    // cout << "delta = " << delta[i] << endl;
        	        step6(i, xt, delta, S);
                }
            } else {
        		if (i > 0) {
                    ksi[i-1] = 0;
        			for (int j = i; j < k; j++)
        				ksi[i-1] += R(i-1, j)*xt[j];
        			dist[i-1] = dist[i] + xidist;
        			i--;
                    step2(i, q, xt, y, delta, ksi, R, S);
        		} else { // lattice point is found (Step 5)
        			radius = dist[0] + xidist;
                    found = true;
        			x = conv_to<vector<int>>::from(xt);
                    i++;
        			step6(i, xt, delta, S);
        		}
        	}
        }
    }
    if (found)
	    return x;
    else {
        log_msg("sphdec: point not found!", "Alert");
        return vector<int>(0);
    }
}

/* Sphere decoder algorithm that considers spherical shaping of the codebook */
vector<int> sphdec_spherical_shaping(const vec &y, const mat &HR, const mat &R, const vector<int> &S,
                                     int &counter, double P, double radius){

    // Step 1
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1;

    // spacing_x = abs(S[0]-S[1]);

    vector<int> x(k); // found point
    vec xt(k), ksi(k), delta(k), dist(k), curr(k), ener(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros(); curr.zeros(); ener.zeros();
    counter = 0;
    
    double xidist = 0.0, xiener = 0.0;
    bool found = false;

    if (!check(radius, HR, y))
        return x;

    step2(i, q, xt, y, delta, ksi, HR, S);

    while (true) {

        if (exit_flag) break; // terminate simulations

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
                step6(i, xt, delta, S);
            }
        } else { // we are inside the sphere

            //xiener = pow((2*xt[i] - q + 1)*R(i,i) + curr[i], 2);
            xiener = pow(xt[i]*R(i,i) + curr[i], 2);

            // we are outside the signal set boundaries
            if (xt[i] < S[0] || xt[i] > S[q - 1] || ener[i] + xiener > P + 10e-6) {
                if ((xt[i] < S[0] && (xt[i] + delta[i]) > S[q - 1]) || (xt[i] > S[q - 1] && (xt[i] + delta[i]) < S[0])){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        step6(i, xt, delta, S);
                    }
                } else 
                    step6(i, xt, delta, S);
            } else {
                if (i > 0) {
                    ksi[i-1] = 0;
                    curr[i-1] = 0;
                    for (int j = i; j < k; j++){
                        ksi[i-1] += HR(i-1, j)*xt[j];
                        // curr[i-1] += (2*xt[j] - q + 1)*R(i-1,j);
                        curr[i-1] += xt[j]*R(i-1,j);
                    }
                    dist[i-1] = dist[i] + xidist;
                    ener[i-1] = ener[i] + xiener;
                    i--;
                    step2(i, q, xt, y, delta, ksi, HR, S);
                } else { // lattice point is found (Step 5)
                    radius = dist[0] + xidist;
                    found = true;
                    x = conv_to<vector<int>>::from(xt);
                    i++;
                    step6(i, xt, delta, S);
                }
            }
        }
    }
    if (found)
        return x;
    else {
        log_msg("sphdec: point not found!", "Alert");
        return vector<int>(0);
    }
}

/* Wrapper function for the sphere decoder to handle complex to real matrix conversion, 
   QR-decomposition and other mappings */
vector<int> sphdec_wrapper(const vector<cx_mat> &bases, const mat Rorig, const cx_mat &H, 
                           const cx_mat &X, const cx_mat &N, const vector<int> &symbset, int &visited_nodes, double radius){

    int n = params["no_of_receiver_antennas"];
    int m = params["no_of_transmit_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    // int q = params["x-PAM"];
    double P = dparams["spherical_shaping_max_power"];

    // cout << "sphdec_wrapper called!" << endl;

    vector<int> x(k);
    cx_mat Y(n, t); //, Ynorm;                 // Helper complex matrices
    mat B(2*t*n, k), G(2*t*m,k), Q, R;         // real matrices
    // vec y(2*t*n), y2;                                // helper and sphdec input vector
    // cx_mat Y;                                  // Helper complex matrices
    // mat B, G, Q, R, Rorig;                     // real matrices
    vec y, y2;                                 // helper and sphdec input vector
    
    Y = H*X + N;                               // Simulated code block that we would receive from MIMO-channel
    // cout << "check 0" << endl;                 
    // Ynorm = (Y + H*basis_sum*(q - 1))*0.5;  // normalize received matrix for the sphere decoder
    y = to_real_vector(Y);                     // convert Y to real vector
    // cout << "check 1" << endl;
    for(int i = 0; i < k; i++){
        B.col(i) = to_real_vector(H*bases[i]); // B = (HX1 HX2 ... HXk) = generator matrix of faded lattice
        // if (P > 0) 
        //     G.col(i) = to_real_vector(bases[i]);   // G = (X1 X2 ... Xk) = generator matrix of original lattice
    }

    // cout << "check 2" << endl;

    // if (P > 0) {
    //     qr_econ(Q, Rorig, G);  // QR-decomposition of G (omits zero rows in Rorig)
    //     process_qr(Q, Rorig);  // Make sure Rorig has positive diagonal elements
    //     Q.zeros();
    // }
    // cout << "check 3" << endl;
    qr_econ(Q, R, B);  // QR-decomposition of B (omits zero rows in R)
    process_qr(Q, R);  // Make sure R has positive diagonal elements

    y2 = Q.st()*y;     // Map y to same basis as R
    // cout << "check 4" << endl;

    if (P < 0)
        x = sphdec(y2, R, symbset, visited_nodes, radius); // Call the actual sphere decoder algorithm
    else
        x = sphdec_spherical_shaping(y2, R, Rorig, symbset, visited_nodes, P, radius);

    // if (x.size() == 0) // point not found
    //     return vector<int>(0);

    // for (int j = 0; j < k; j++)
    //     x[j] = 2*x[j] - q + 1;

    return x;
}