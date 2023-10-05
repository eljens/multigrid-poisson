#ifndef INPUT_PARSER
#define INPUT_PARSER

#include <iostream>
#include <string>
#include "definitions.h"
#include "settings.h"
#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::endl;
using std::string;
using std::stoi;
using std::stod;
using std::ostream;

namespace Poisson{
    Settings parser(int argc, char * argv[]);
}

#endif
