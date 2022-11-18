#include "include/parser.h"
#include "include/domain.h"
#include "include/definitions.h"
#include "include/problem_definition.h"
#include "include/jacobi.h"
#include <iostream>
#include <cstdlib>
#include <cassert>

using std::cout;
using std::endl;

int main(int argc, char * argv[]){
    Settings settings = parser(argc,argv);
    Domain<double_t> domain(settings,false);

    cout << "Created domain" << endl;
    const double_t omega = 2.0/3.0;

    domain.init(&ufun,&ffun,&dudxfun,&dudyfun);

    cout << "Initialized domain" << endl;

    domain.to_device();

    cout << "To device finished" << endl;

    for(int_t i=0;i<settings.maxiter;i++){
        jacobi<double_t>(domain,omega);
    }

    domain.to_host();

    domain.save();

    domain.save_halo();

    double_t err = 0.0;

    #pragma omp parallel for collapse(3) reduction(max:err)
    for (int_t i = 0;i<domain.u->shape[0];i++){
        for (int_t j = 0;j<domain.u->shape[1];j++){
            for (int_t k = 0;k<domain.u->shape[2];k++){
                double_t tmp = abs(domain.u->at[domain.u->idx(i,j,k)]
                    -ufun(settings.origin[0]+((double_t) i)*settings.h,
                    settings.origin[1]+((double_t) j)*settings.h,
                    settings.origin[2]+((double_t) k)*settings.h));
                err = std::max(err,tmp);
            }
        }
    }
 
    cout << "Maximal error: " << err << endl;
    return EXIT_SUCCESS;
}