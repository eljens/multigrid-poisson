#ifndef INPUT_SETTINGS
#define INPUT_SETTINGS

#include <iostream>
#include <string>
#include <omp.h>
#include "definitions.h"

using std::cout;
using std::endl;
using std::ostream;

namespace Poisson{
    class Settings {
        public:
            uint_t levels = 1;
            uint_t dims[3] = {3,3,3};
            uint_t maxiter = 100;
            double_t tolerance = 1e-6;
            int_t host = omp_get_initial_device();
            int_t dev = omp_get_default_device();
            double_t lengthx = 3.1415; // Length in x dimension
            double_t h = 0.33; 
            double_t origin[3] = {-0.4, 0.21, 1.23};
            bool print_result = false;
            bool write_final_stats = false;
            std::string stats_file = "";

            Settings();

            Settings(Settings & settings);

            ~Settings();
    };

    ostream& operator<<(ostream& os, const Settings& settings);
}
#endif
