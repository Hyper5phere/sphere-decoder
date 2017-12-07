/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : sphdec.cpp                                                                            *
 * Project     : Schnorr-Euchnerr sphere decoder simulation for space-time lattice codes               *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++11                                                                                 *
 * Description : The core algorithm implementations for the simulation + wrapper function for them     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ARMA_NO_DEBUG /* disable Armadillo bound checks for addiotional speed */

#include <iostream>
#include <armadillo> /* linear algebra library */
#include <vector>
#include <complex>
#include <string>

#include "sphdec.hpp"
#include "misc.hpp"
#include "algorithms.hpp"

using namespace std;
using namespace arma;


/* Decision feedback equalization on xt[i] 
   i.e. Use the Babai nearest plane algorithm as a starting point in dimension i */
inline void dfe(int i, int q, vec &xt, const vec &y, vec &delta, const vec &ksi, const mat &R, const vector<int> &S) {
    /* calculate the "Babai coeffient" and "round" it to nearest symbol in the symbol set */
    xt[i] = nearest_symbol((y[i]-ksi[i])/R(i,i), S);
    if (xt[i] <= S[0]) { /* only feasible way is positive direction */
        delta[i] = 2;
    } else if (xt[i] >= S[q - 1]) { /* only feasible way is negative direction */
        delta[i] = -2;
    } else { /* there's room to enumerate within signal set */
        delta[i] = 2*sesd_sign(y[i]-ksi[i]-R(i,i)*xt[i]);
    }
}

/* Schnorr-Euchner enumeration step for the sphere decoder i.e. zig-zag around the dfe point.
   Makes sure we loop over the points in dimension i in non descending order.
   This enables early stopping criteria when we find a point in dimension i outside the radius, 
   i.e. no points in that sub search tree are inside the sphere because they can only be further away or at the same distance. */
inline void se_enum(int i, vec &xt, vec &delta) {
    /* add delta to xt and pick new delta (difference to next symbol in S we want to try) */
    xt[i] = xt[i] + delta[i];
    delta[i] = -delta[i] - 2*sesd_sign(delta[i]);
}

/* Perform some basic checks to ensure that the sphere decoder works as intended */
inline bool check(double radius, const mat &R, const vec &y) {
    if (radius <= 0){
        log_msg("sphdec: negative squared initial radius given!", "Error");
        log_msg("aborting simulation round...");
        return false;
    }

    if (R.n_rows != R.n_cols) {
        log_msg("sphdec: R is not a square matrix!", "Error");
        log_msg("aborting simulation round...");
        return false;
    }

    if (y.n_elem != R.n_cols) {
        log_msg("sphdec: vector y dimension mismatch!", "Error");
        log_msg("aborting simulation round...");
        return false;
    }
    return true;
}


/* Basic sphere decoder algorithm */
vector<int> sphdec(const vec &y, const mat &R, const vector<int> &S, int &counter, double radius){

    /* Initialize */
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1; /* We start from dimension k-1 and iterate all the way down to dimension 0 */

    vector<int> x(k); /* point to decode */
    vec xt(k), ksi(k), delta(k), dist(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros();
    counter = 0; /* counts how many search tree nodes (loop iterations) we went through */
    
    double xidist = 0.0;
    bool found = false;

    if (!check(radius, R, y))
        return vector<int>(0);

    /* Initialize xt[k-1] and delta[k-1] */
    dfe(i, q, xt, y, delta, ksi, R, S);

    while (!exit_flag) {

        counter++;
        /* Step 3. */
        xidist = pow(y[i]-ksi[i]-R(i,i)*xt[i], 2);
        
        /***** Uncomment line below to debug *****/
        // cout << i << ": " << vec2str(xt, k) << ", xidist = " << xidist + dist[i] << ", C = " << radius << endl;

        if (radius < dist[i] + xidist) { // current point xt is outside the sphere
            // Step 4.
            if (i == k-1) {
                break;
            } else {
                // Step 6.
                i++;
                se_enum(i, xt, delta);
            }
        } else { // we are inside the sphere
            if (xt[i] < S[0] || xt[i] > S[q - 1]){ // we are outside the signal set boundaries
                if ((xt[i] < S[0] && (xt[i] + delta[i]) > S[q - 1]) || (xt[i] > S[q - 1] && (xt[i] + delta[i]) < S[0])){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        se_enum(i, xt, delta);
                    }
                } else {
                    se_enum(i, xt, delta);
                }
            } else {
                if (i > 0) {
                    ksi[i-1] = 0;
                    for (int j = i; j < k; j++)
                        ksi[i-1] += R(i-1, j)*xt[j];
                    dist[i-1] = dist[i] + xidist;
                    i--;
                    dfe(i, q, xt, y, delta, ksi, R, S);
                } else { // lattice point is found (Step 5)
                    radius = dist[0] + xidist;
                    found = true;
                    x = conv_to<vector<int>>::from(xt);
                    i++;
                    se_enum(i, xt, delta);
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

/* Sphere decoder algorithm that considers the spherical shaping of the codebook */
vector<int> sphdec_spherical_shaping(const vec &y, const mat &HR, const mat &R, const vector<int> &S,
                                     int &counter, double P, double radius){

    /* Initialize */
    int k = params["no_of_matrices"];
    int q = params["x-PAM"];
    int i = k-1; /* We start from dimension k-1 and iterate all the way down to dimension 0 */

    vector<int> x(k); /* point to decode */
    vec xt(k), ksi(k), delta(k), dist(k), curr(k), ener(k);
    xt.zeros(); ksi.zeros(); delta.zeros(); dist.zeros(); curr.zeros(); ener.zeros();
    counter = 0; /* counts how many search tree nodes (loop iterations) we went through */
    
    double xidist = 0.0, xiener = 0.0;
    bool found = false;

    if (!check(radius, HR, y))
        return vector<int>(0);

    /* Initialize xt[k-1] and delta[k-1] */
    dfe(i, q, xt, y, delta, ksi, HR, S);

    while (!exit_flag) {

        counter++;
        // Step 3.
        xidist = pow(y[i]-ksi[i]-HR(i,i)*xt[i], 2);
        
        /***** Uncomment line below to debug *****/
        // cout << i << ": " << vec2str(xt, k) << ", xidist = " << xidist + dist[i] << ", C = " << radius << endl;

        if (radius < dist[i] + xidist) { /* current point xt is outside the sphere */
            // Step 4.
            if (i == k-1) {
                break;
            } else {
                // Step 6.
                i++;
                se_enum(i, xt, delta);
            }
        } else { /* we are inside the sphere */

            xiener = pow(xt[i]*R(i,i) + curr[i], 2);

            // we are outside the signal set boundaries
            if (xt[i] < S[0] || xt[i] > S[q - 1] || ener[i] + xiener > P + 10e-6) {
                if ((xt[i] < S[0] && (xt[i] + delta[i]) > S[q - 1]) || (xt[i] > S[q - 1] && (xt[i] + delta[i]) < S[0])){
                    // Step 4.
                    if (i == k-1) {
                        break;
                    } else {
                        i++;
                        se_enum(i, xt, delta);
                    }
                } else 
                    se_enum(i, xt, delta);
            } else {
                if (i > 0) {
                    ksi[i-1] = 0;
                    curr[i-1] = 0;
                    for (int j = i; j < k; j++){
                        ksi[i-1] += HR(i-1, j)*xt[j];
                        curr[i-1] += xt[j]*R(i-1,j);
                    }
                    dist[i-1] = dist[i] + xidist;
                    ener[i-1] = ener[i] + xiener;
                    i--;
                    dfe(i, q, xt, y, delta, ksi, HR, S);
                } else { // lattice point is found (Step 5)
                    radius = dist[0] + xidist;
                    found = true;
                    x = conv_to<vector<int>>::from(xt);
                    i++;
                    se_enum(i, xt, delta);
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
                           const cx_mat &X, const cx_mat &N, const vector<int> &symbset, int &visited_nodes, double radius) {

    /* read simulation parameters */
    int n = params["no_of_receiver_antennas"];
    int t = params["time_slots"];
    int k = params["no_of_matrices"];
    double P = dparams["spherical_shaping_max_power"];

    vector<int> x(k);                          /* decoded point */
    cx_mat Y(n, t);                            /* received code block */
    mat B(2*t*n, k), Q, R;                     /* real matrices for QR-decompostion */
    vec y(2*t*n);                              /* raw input vector for the sphere decoder (unmapped) */
    vec y2;                                    /* input vector for the sphere decoder (mapped to same space a R) */
    
    Y = H*X + N;                               /* Calculate simulated code block that we would receive */          
    y = to_real_vector(Y);                     /* convert Y to real vector */

    for (int i = 0; i < k; i++)
        B.col(i) = to_real_vector(H*bases[i]); /* B = (HX1 HX2 ... HXk) = generator matrix of faded lattice */

    qr_econ(Q, R, B);                          /* QR-decomposition of B (omits zero rows in R) */
    process_qr(Q, R);                          /* Make sure R has positive diagonal elements */

    y2 = Q.st()*y;                             /* Map y to same basis as R */

    /* decide which sphere decoder algorithm to use */
    if (P <= 0)
        x = sphdec(y2, R, symbset, visited_nodes, radius); 
    else
        x = sphdec_spherical_shaping(y2, R, Rorig, symbset, visited_nodes, P, radius);

    return x;
}