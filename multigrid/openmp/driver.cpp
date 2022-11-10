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

    /*
    DeviceArray<double_t>A(omp_get_default_device(),{10,3,2});
    for (int i=0;i<10;i++){
        for (int j=0;j<3;j++){
            for (int k=0;k<2;k++){
                A.at[A.idx({i,j,k})] = i*3*2 + j*2 + k;
            }
        }
    }
    A.to_device();
    double_t * Ad = A.devptr;
    #pragma omp target teams distribute parallel for schedule(static) is_device_ptr(Ad) device(A.device)
    for (int i=0;i<10;i++){
        for (int j=0;j<3;j++){
            for (int k=0;k<2;k++){
                Ad[A.idx({i,j,k})] = 2*Ad[A.idx({i,j,k})];
            }
        }
    }
    A.to_host();
    for (int i=0;i<A.size;i++){
        assert(A.at[i]==2*i);
    }
    cout << "DeviceArray works as intended" << endl;
    */
    cout << "Yas" << endl;
    return EXIT_SUCCESS;
}