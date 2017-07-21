#include <iostream>
#include <mutex>
#include <ctime>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

#include "misc.hpp"
#ifdef PLOTTING
#include "gnuplot-iostream.hpp"
#endif

using namespace std;

/* object for parallel computing synchronazation */
std::mutex log_mutex;
std::mutex output_mutex;

/* Generates a standard datetime string of current local time */
string time_str(){
    struct tm *timeinfo;
    time_t t = time(nullptr);
    timeinfo = localtime(&t);
    char timestamp[50]; 
    strftime(timestamp, 50, "%d-%m-%Y %T", timeinfo);
    return string(timestamp);
}

/* function used for logging, should be thread safe */
void log_msg(const string msg, const string lvl){
        string prefix = time_str();
        prefix += " | [" + lvl + "]\t";
        lock_guard<mutex> lock(log_mutex); // make sure other threads don't write to log file simultaneosly
        ofstream logfile(filenames["log"], ios_base::app);
        if (msg.compare("-start-") == 0)
            logfile << endl;
        else {
            logfile << prefix << msg << endl;
            cout << prefix << msg << endl;
        }
        logfile.close();
}

/* Makes the user input somewhat more readable */
void clean_input(string &input){
    // remove comments
    string::size_type start = input.find("//", 0);
    if (start != string::npos)
        input = input.substr(0, start);
    // remove white spaces
    input.erase(remove_if(input.begin(), input.end(), ::isspace), input.end());
}

/* Generates the output filename from current time and the used bases file name */
void create_output_filename(){
    string bf = filenames["bases"];
    size_t a = bf.find("/")+1;
    size_t b = bf.find(".");
    string name = bf.substr(a, b-a);
    filenames["output"] = string("output/") + time_str() + string(" ") + name + string(" output.csv");
    // return string("output/") + time_str() + string(" ") + name + string(" output.csv");
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
    #ifdef PLOTTING
    Gnuplot gp;
    gp << "set terminal x11\n"; // Open plots in new GUI window
    if (logscale)
        gp << "set logscale y 10\n";  // Set y axis to logscale of base 10 (SNR is already in log10 scale)
    gp << "set datafile separator \",\"\n";
    gp << "set xlabel \"" << xlabel << "\"\n";
    gp << "set ylabel \"" << ylabel << "\"\n";
    gp << "set grid\n";
    gp << "plot '" << filenames["output"] << "' using " << xcol << ":" << ycol << " with linespoints ls 3\n";
    #endif
}

void output_data(parallel_vector<string> &output){
    if (output.size() > 1){
        sort(output.begin()+1, output.end(), snr_ordering); // sort vector by SNR (ascending order)
        log_msg("Printing simulation output to '" + filenames["output"] + "'...");
        output_csv(output);
        #ifdef PLOTTING // Draw plots with Gnuplot if plotting is enabled
        if (params["plot_results"] > 0){
            log_msg("Drawing plots...");
            plot_csv(1, 4, "SNR (dB)", "BLER (%)", true);
            plot_csv(1, 5, "SNR (dB)", "Average Complexity (# visited points)", false);
        }
        #endif
    }
}
