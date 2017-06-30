#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <map>
#include <vector>
#include <mutex>
#include <set>

/* default filenames */
extern std::string options_filename;
extern std::string basis_filename;
extern std::string log_filename;

/* storage for simulation parameters */
extern std::map<std::string, int> params;

/* thread safe vector data structure */
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

/* thread safe set data structure */
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
void log_msg(const std::string msg = "-exit-", const std::string lvl = "Info");
void clean_input(std::string &input);
std::string create_output_filename(void);
void output_csv(const std::string &filename, const parallel_vector<std::string> &line);
bool snr_ordering(std::string &a, std::string &b);

/* Makes a string representation out of any basic vector type */
template <typename T>
inline std::string vec2str(T vec, size_t size){
	std::string str = "{";
	for (size_t i = 0; i < (size-1); i++)
		str += std::to_string(vec[i]) + ", ";
	return str + std::to_string(vec[size-1]) + "}";
}

#endif /* MISC_HPP */