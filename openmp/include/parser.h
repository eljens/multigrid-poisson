#ifndef INPUT_PARSER
#define INPUT_PARSER

#include <iostream>
#include <string>
#include <omp.h>
#include "definitions.h"
#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::endl;
using std::string;
using std::stoi;
using std::stod;
using std::ostream;

class Settings {
    public:
    uint_t dims[3] = {0,0,0};
    uint_t maxiter = 100;
    double_t tolerance = 1e-6;
    int_t host = omp_get_initial_device();
    int_t dev = omp_get_default_device();
    uint_t numdev = omp_get_num_devices();
    double_t h = 0.01; 
    double_t origin[3] = {0.5, -0.5, -1.25};
};

Settings parser(int argc, char * argv[]);

ostream& operator<<(ostream& os, const Settings& settings);

#endif