#ifndef MISC_HPP
#define MISC_HPP

#define _GLIBCXX_USE_CXX11_ABI 0 // fixes some string related errors
#define ARMA_NO_DEBUG // for speed
// #define PLOTTING // enable plotting (requires boost C++ library)

#include <string>
#include <map>
#include <vector>
#include <mutex>
#include <set>
#include <armadillo>
#include <iomanip>
#include <sstream>

/* storage for filenames */
extern std::map<std::string, std::string> filenames;

/* storage for simulation integer parameters */
extern std::map<std::string, int> params;

/* storage for simulation double parameters */
extern std::map<std::string, double> dparams;

/* storage for simulation string parameters */
extern std::map<std::string, std::string> sparams;

/* flag which indicates early stop of the simulation */
extern bool exit_flag;

/* thread safe std::vector data structure */
template <typename T>
class parallel_vector : public std::vector<T> {
	public:
		parallel_vector() : std::vector<T>(){};
		void append(T item)	{
		    std::lock_guard<std::mutex> lock(this->m_);
		    this->push_back(item);
		};
	private:
		std::mutex m_;
};

/* thread safe std::set data structure */
template <typename T>
class parallel_set : public std::set<T> {
	public:
		parallel_set() : std::set<T>(){};
		void par_insert(T item)	{
		    std::lock_guard<std::mutex> lock(this->m_);
		    this->insert(item);
		};
	private:
		std::mutex m_;
};

/* function prototypes */
std::string time_str(void);
void log_msg(const std::string msg = "-start-", const std::string lvl = "Info");
void clean_input(std::string &input);
void create_output_filename(void);
void output_csv(const parallel_vector<std::string> &line);
bool snr_ordering(std::string &a, std::string &b);
void plot_csv(int xcol, int ycol, const std::string &xlabel, const std::string &ylabel, bool logscale);
void output_data(parallel_vector<std::string> &output);
void output_real_matrix(const std::string &filepath, const arma::mat &A);
void output_complex_matrix(const std::string &filepath, const arma::cx_mat &A);

/* Makes a string representation out of any basic vector type */
template <typename T>
inline std::string vec2str(T vec, size_t size){
	std::string str = "{";
	for (size_t i = 0; i < (size-1); i++)
		str += std::to_string(vec[i]) + ", ";
	return str + std::to_string(vec[size-1]) + "}";
}

/* converts a floating point value to a string with custom precision
 * credits go to: https://stackoverflow.com/a/16606128 
 */
template <typename T>
inline std::string float2str(const T value, int precision = 6){
	std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

/* returns a string representation of an Armadillo matrix object */
template <typename T>
inline std::string mat2str(const T A){
	std::ostringstream out;
    out << std::endl << A;
    return out.str();
}

#endif /* MISC_HPP */