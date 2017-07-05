#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

void create_config(void);
void configure(void);
std::vector<arma::cx_mat> read_matrices(void);

#endif /* CONFIG_HPP */