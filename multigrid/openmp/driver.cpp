#import "array.h"
#import "parser.h"
#import "domain.h"
#include "definitions.h"
#include "problem_definition.h"
#import <iostream>
#import <cstdlib>
#import <cassert>

using std::cout;
using std::endl;

int main(int argc, char * argv[]){
    Settings settings = parser(argc,argv);
    Domain<double_t> domain(settings);

    domain.init(&ufun,&ffun,&dudxfun,&dudyfun);

    domain.to_device();

    domain.to_host();

    domain.save_result();

    
 
    cout << "Yas" << endl;
    return EXIT_SUCCESS;
}