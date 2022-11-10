#ifndef INPUT_PARSER
#define INPUT_PARSER

#include <iostream>
#include <string>
#include <omp.h>
#include "definitions.h"

using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::stoi;
using std::stod;

class Settings {
    public:
    uint_t dims[3] = {0,0,0};
    uint_t maxiter = 100;
    double_t tolerance = 1e-6;
    int_t host = omp_get_initial_device();
    uint_t numdev = omp_get_num_devices();
};

Settings parser(int argc, char * argv[]);

std::ostream& operator<<(std::ostream& os, const Settings& settings);

#endif