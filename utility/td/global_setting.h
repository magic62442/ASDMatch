#ifndef GLOBAL_SETTING_HZY_H
#define GLOBAL_SETTING_HZY_H

// #define MEMORIZATION
#define TIMEOUT
#define PREPROCESS
#define PRINTSTATENUMBER
// #define DIVBICCP
//#define IDV_TEST
// #define ENUM
// #define CHECKF
// #define DVO
#define DPOD
#define STATISTIC
#define OUTPUTTD
// #define SLMPV

// #define DIVAPPROX
// #define BYCOVER
// #define CALC_WIDTH
// #define MULTITHREAD
#define MULTIPROGRESS
#define HLOG

// #define LOW_DELTA2
#define LOW_GAMMAR
//#define LOW_MMDP
#define UP_MINW

#include <time.h>
#include <string.h>

const double eps = 1e-2; 
const double inf = 1e7;
const size_t tau = 5; // n <= tau using enumerate method
const size_t lambda = 20; // n >= lamda using approximate method
const size_t mem_lim = 20;
const size_t UINT64NUM = 20;
const size_t MAXBITNUMBER = UINT64NUM * 64;
const size_t MAXELENUM = MAXBITNUMBER;
const size_t TIMEOUTSEC = 1;
const clock_t time_out = CLOCKS_PER_SEC * TIMEOUTSEC;


const char * const result_prefix = "../result/";
#ifdef UP_MINW
const char * const ub = "minw_";
#else 
const char * const ub = "noub_";
#endif

#ifdef LOW_DELTA2
const char * const lb = "delta2_";
#else 
#ifdef LOW_GAMMAR
const char * const lb = "gammaR_";
#else
#ifdef LOW_MMDP
const char * const lb = "mmd+_";
#else
const char * const lb = "nolb_";
#endif
#endif
#endif

// std::string const & temp = std::string(result_prefix) + std::string(ub) + std::string(lb) + std::to_string(TIMEOUTSEC);
// const char * const res_dir = (std::string(result_prefix) + std::string(ub) + std::string(lb) + std::to_string(TIMEOUTSEC)).c_str();
// #define res_dir (std::string(result_prefix) + std::string(ub) + std::string(lb) + std::to_string(TIMEOUTSEC)).c_str()
const char * const data_dir = "../graph";
const char * const idv_file = "./idv.txt";


const size_t worker_num = 4;

#endif
