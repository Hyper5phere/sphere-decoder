#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <map>

/* default filenames */
extern std::string options_filename;
extern std::string basis_filename;
extern std::string output_filename;
extern std::string log_filename;

/* storage for simulation parameters */
extern std::map<std::string, int> params;

std::string time_str(void);
void log_msg(const std::string msg = "-exit-", const std::string lvl = "Info");
void clean_input(std::string &input);

#endif /* MISC_HPP */