/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Filename    : config.hpp                                                                            *
 * Project     : Schnorr-Euchnerr sphere decoder simulation for space-time lattice codes               *
 * Authors     : Pasi Pyrrö, Oliver Gnilke                                                             *
 * Version     : 1.0                                                                                   *
 * Copyright   : Aalto University ~ School of Science ~ Department of Mathematics and Systems Analysis *
 * Date        : 17.8.2017                                                                             *
 * Language    : C++11                                                                                 *
 * Description : Function prototypes for config.cpp                                                    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

void create_config(void);
void configure(void);
std::vector<arma::cx_mat> read_matrices(const std::string &filepath);
std::unordered_map<int, int> read_error_requirements(const std::string &filepath);

#endif /* CONFIG_HPP */