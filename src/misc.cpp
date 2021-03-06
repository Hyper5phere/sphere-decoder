/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : misc.cpp                                                                              *
 * Project     : Planewalker - Schnorr-Euchner sphere decoder simulation for space-time lattice codes  *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++ (2011 or newer standard)                                                          *
 * Description : contain miscellaneous utilities for the sphere decoder program, mainly I/O stuff      *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ARMA_NO_DEBUG /* disable Armadillo bound checks for addiotional speed */

#include <iostream>
#include <mutex>
#include <ctime>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <armadillo>

#include "misc.hpp"
#ifdef PLOTTING
#include "gnuplot-iostream.hpp" /* Stand alone header file for Gnuplot API (requires Boost library though) */
#endif

using namespace std;
using namespace arma;

/* objects for parallel computing synchronazation */
std::mutex log_mutex;
std::mutex output_mutex;

string red = "\e[1;31m";
string green = "\e[1;32m";
string yellow = "\e[1;33m";
string def = "\e[0m";

/* Generates a standard datetime string of current local time */
string time_str(){
    struct tm *timeinfo;
    time_t t = time(nullptr);
    timeinfo = localtime(&t);
    char timestamp[50]; 
    strftime(timestamp, 50, "%d-%m-%Y %T", timeinfo);
    return string(timestamp);
}

/* function used for logging, should be thread safe, outputs to terminal and into a txt file
 *
 * Defined values for lvl argument:
 * "Raw"   : msg is printed as is without prefix or color
 * "Info"  : default value, has prefix but no color, used for basic logging
 * "Alert" : Something unexpected happened, uses yellow color
 * "Error" : An error occurred and the program needs to terminate, uses red color
 */
void log_msg(const string msg, const string lvl) {
    string prefix = "";
    string color = "";
    if (lvl.compare("Raw") != 0) { /* raw input contains no prefix */
        prefix += time_str();
        prefix += " | [" + lvl + "]\t";
    }
    if (lvl.compare("Alert") == 0){
        color = yellow;
    } 
    if (lvl.compare("Error") == 0){
        color = red;
    }
    lock_guard<mutex> lock(log_mutex); /* make sure other threads don't write to log file simultaneosly */
    ofstream logfile(filenames["log"], ios_base::app);
    if (msg.compare("-start-") == 0)   /* just writes newline to log file to indicate new program run */
        logfile << endl;
    else {
        logfile << prefix << msg << endl; /* log txt file */
        cout << color << prefix << msg << def << endl;    /* stdout */
    }
    logfile.close();
}

/* Makes the user input somewhat more readable */
void clean_input(string &input){
    /* remove comments from string */
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    /* remove white spaces from string */
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
}

/* removes all characters c from input string */
void remove_from_string(string &input, const char c){
    input.erase(remove(input.begin(), input.end(), c), input.end());
}

/* Generates the output filename from current time and the used bases file name */
void create_output_filename(){
    string bf = filenames["bases"];
    size_t a = bf.find("/")+1;
    size_t b = bf.find(".");
    string name = bf.substr(a, b-a);
    filenames["output"] = string("output/") + time_str() + string(" ") + name + string(" output.csv");
}

/* Writes a vector of strings (lines) in a csv file */
void output_csv(const parallel_vector<string> &lines){
    ofstream csv(filenames["output"]);
    for (const auto &line : lines){
        csv << line << endl;
    }
    csv.close();
}

/* Used to compare the SNR values of each output line */
bool snr_ordering(string &a, string &b){
    /* Example: 
     * a = "16,15.976382,10000,9.760000"
     * b = "18,17.976382,10000,10.160000"
     * --> i = 16, j = 18
     * --> a comes before b
     */
    int i = strtol(a.substr(0, a.find(",")).c_str(), NULL, 10);
    int j = strtol(b.substr(0, b.find(",")).c_str(), NULL, 10);
    return (i < j);
}


/* Uses gnuplot stream to plot the columns of the output csv file */
void plot_csv(int xcol, int ycol, const string &xlabel, const string &ylabel, bool logscale){
    #ifdef PLOTTING /* This function does nothing if PLOTTING is not defined when compiled */
    Gnuplot gp;
    gp << "set terminal x11\n"; /* Open plots in a new GUI window by default, can also be saved in png file for example */
    if (logscale)
        gp << "set logscale y 10\n";  /* Set y axis to logscale of base 10 (SNR is already in log10 scale) */
    gp << "set datafile separator \",\"\n";
    gp << "set xlabel \"" << xlabel << "\"\n";
    gp << "set ylabel \"" << ylabel << "\"\n";
    gp << "set grid\n";
    /* plot using the data from output csv file */
    gp << "plot '" << filenames["output"] << "' using " << xcol << ":" << ycol << " with linespoints ls 3\n";
    #endif
}

/* Prints simulation results in a file and optionally plots the results */
void output_data(parallel_vector<string> &output){
    if (output.size() > 1){
        sort(output.begin()+1, output.end(), snr_ordering); /* sort vector by SNR (ascending order) */
        log_msg("Printing simulation output to '" + filenames["output"] + "'...");
        output_csv(output);
        #ifdef PLOTTING /* Draw plots with Gnuplot if plotting is enabled */
        if (params["plot_results"] > 0){
            log_msg("Drawing plots...");
            plot_csv(2, 7, "SNR (dB)", "BLER (%)", true);
            plot_csv(2, 3, "SNR (dB)", "Average Complexity (# visited points)", false);
        }
        #endif
    }
}


/* outputs real matrix into a file in mathematica format */
void output_real_matrix(const string &filepath, const mat &A, bool append){
    ofstream mfile;
    /* selects write mode: either overwrite the old file or append to it */
    if (!append)
        mfile = ofstream(filepath);
    else
        mfile = ofstream(filepath, ios_base::app);
    mfile << "{";
    for (auto i = 0u; i < A.n_rows; i++){
        mfile << "{";
        for (auto j = 0u; j < A.n_cols; j++){
            mfile << float2str(A(i, j), 16);
            if (j < A.n_cols - 1)
                mfile << ", ";
        }
        mfile << "}";
        if (i < A.n_rows - 1)
            mfile << "," << endl;
    }
    mfile << "}";
    if (append) mfile << endl << endl;
    mfile.close();
}


/* outputs complex matrix into a file in mathematica format */
void output_complex_matrix(const string &filepath, const cx_mat &A, bool append){
    ofstream mfile;
    string real, imag;
    /* selects write mode: either overwrite the old file or append to it */
    if (!append)
        mfile = ofstream(filepath);
    else
        mfile = ofstream(filepath, ios_base::app);
    mfile << "{";
    for (auto i = 0u; i < A.n_rows; i++){
        mfile << "{";
        for (auto j = 0u; j < A.n_cols; j++){
            real = float2str(A(i, j).real(), 16);
            imag = float2str(A(i, j).imag(), 16);
            if (imag[0] == '-')
                mfile << real << " - " << imag.substr(1) << " *I";
            else
                mfile << real << " + " << imag << " *I";
            if (j < A.n_cols - 1)
                mfile << ", ";
        }
        mfile << "}";
        if (i < A.n_rows - 1)
            mfile << "," << endl;
    }
    mfile << "}";
    if (append) mfile << endl << endl;
    mfile.close();
}
