#ifndef MISC_HPP
#define MISC_HPP

#define _GLIBCXX_USE_CXX11_ABI 0 // fixes some string related errors
// #define ARMA_NO_DEBUG // for speed
#define PLOTTING // enable plotting (requires boost C++ library)

#include <string>
#include <map>
#include <vector>
#include <mutex>
#include <set>

/* storage for filenames */
extern std::map<std::string, std::string> filenames;

/* storage for simulation integer parameters */
extern std::map<std::string, int> params;

/* storage for simulation double parameters */
extern std::map<std::string, double> dparams;

/* storage for simulation string parameters */
extern std::map<std::string, std::string> sparams;

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

/* Makes a string representation out of any basic vector type */
template <typename T>
inline std::string vec2str(T vec, size_t size){
	std::string str = "{";
	for (size_t i = 0; i < (size-1); i++)
		str += std::to_string(vec[i]) + ", ";
	return str + std::to_string(vec[size-1]) + "}";
}

#endif /* MISC_HPP */