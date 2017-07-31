#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

void create_config(void);
void configure(void);
std::vector<arma::cx_mat> read_matrices(const std::string &filepath);
std::map<int, int> read_error_requirements(const std::string &filepath);

#endif /* CONFIG_HPP */