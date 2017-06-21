#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>
#include <mutex>
#include <ctime>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>

/* default filenames */
extern std::string options_filename;
extern std::string basis_filename;
extern std::string output_filename;
extern std::string log_filename;

/* storage for simulation parameters */
extern std::map<std::string, int> params;

// std::mutex log_mutex;

namespace misc {
	void log_msg(const std::string &);
}
void clean_input(std::string &);

#endif /* MISC_HPP */