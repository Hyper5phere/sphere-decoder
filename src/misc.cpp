#include <iostream>
#include <mutex>
#include <ctime>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

#include "misc.hpp"

using namespace std;

/* object for parallel computing synchronazation */
mutex log_mutex;

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
    ofstream logfile(log_filename, ios_base::app);
    if (msg.compare("-exit-") == 0)
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