#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <map>
#include <vector>

/* default filenames */
extern std::string options_filename;
extern std::string basis_filename;
extern std::string output_filename;
extern std::string log_filename;

/* storage for simulation parameters */
extern std::map<std::string, int> params;

/* function prototypes */
std::string time_str(void);
void log_msg(const std::string msg = "-exit-", const std::string lvl = "Info");
void clean_input(std::string &input);
std::string create_output_filename(void);
void output_csv_line(const std::string &filename, const std::vector<std::string> &line);

/* Makes a string representation out of any basic vector type */
template <typename T>
inline std::string vec2str(T vec, size_t size){
	std::string str = "{";
	for (size_t i = 0; i < (size-1); i++)
		str += std::to_string(vec[i]) + ", ";
	return str + std::to_string(vec[size-1]) + "}";
}

#endif /* MISC_HPP */