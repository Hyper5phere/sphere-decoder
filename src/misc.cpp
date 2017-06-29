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
mutex output_mutex;

void parallel_vector::append(string item){
    lock_guard<mutex> lock(this->m_);
    this->push_back(item);
}

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

string create_output_filename(){
    string bf = basis_filename;
    size_t a = bf.find("/")+1;
    size_t b = bf.find(".");
    string name = basis_filename.substr(a, b-a);
    return string("output/") + time_str() + string(" ") + name + string(" output.csv");
}

void output_csv(const string &filename, const parallel_vector &lines){
    ofstream csv(filename);
    for (const auto &line : lines){
        csv << line << endl;
    }
    csv.close();
}

bool snr_ordering(string &a, string &b){
    int i = strtol(a.substr(0, a.find(",")).c_str(), NULL, 10);
    int j = strtol(b.substr(0, b.find(",")).c_str(), NULL, 10);
    return (i < j);
}