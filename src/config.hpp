#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

void create_config(const std::string);
void configure(const std::string);
std::vector<arma::cx_mat> read_matrices();

#endif /* CONFIG_HPP */