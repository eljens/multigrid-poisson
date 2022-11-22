#include "include/parser.h"
#include "include/domainsettings.h"
#include "include/domain.h"
#include "include/definitions.h"
#include "include/problem_definition.h"
#include "include/jacobi.h"
#include "include/residual.h"
#include "include/injection.h"
#include "include/trilinearinterpolation.h"
#include "include/grid.h"
#include "include/vcycle.h"
#include <iostream>
#include <cstdlib>
#include <cassert>

using std::cout;
using std::endl;

int main(int argc, char * argv[]){
    bool is_dirichlet = false;
    Settings settings = parser(argc,argv);
    DomainSettings domainsettings(settings,0);
    Grid grid(settings,settings.levels);
    const double_t omega = 2.0/3.0;

    // Making array of domains
    Domain<double_t> * domains[settings.levels];
    for (uint_t l = 0;l<settings.levels;l++){
        domains[l] = new Domain<double_t>(grid.domainsettings[l],is_dirichlet);
        if (l==0){
            domains[l]->init(&ufun,&ffun,&dudxfun,&dudyfun);
        }
        else {
            domains[l]->init_zero();
        }
        domains[l]->to_device();
    }

    Injection<double_t> injection;
    TrilinearInterpolation<double_t> trilinearinterpolation;
    double_t start = omp_get_wtime();
    for(uint_t i = 0;i<settings.maxiter;i++){
        Vcycle<double_t>(domains,injection,trilinearinterpolation,omega,0,settings.levels);
        if (i % 20 == 19){
            cout << "Finished iteration " << i+1 << endl;
        }
    }
    cout << endl;
    cout << "It took " << omp_get_wtime()-start << " seconds to run " <<settings.maxiter << " Vcycles"  <<endl;

    residual<double_t>(*domains[0]);

    domains[0]->to_host();

    if (settings.print_result){
        domains[0]->save();
    }

    //domains[0]->save_halo();

    //domains[0]->save_restrict_prolong(injection,trilinearinterpolation);

    double_t err = 0.0;

    #pragma omp parallel for collapse(3) reduction(max:err)
    for (int_t i = 0;i<domains[0]->u->shape[0];i++){
        for (int_t j = 0;j<domains[0]->u->shape[1];j++){
            for (int_t k = 0;k<domains[0]->u->shape[2];k++){
                double_t tmp = abs(domains[0]->u->at[domains[0]->u->idx(i,j,k)]
                    -ufun(settings.origin[0]+((double_t) i)*settings.h,
                    settings.origin[1]+((double_t) j)*settings.h,
                    settings.origin[2]+((double_t) k)*settings.h));
                err = std::max(err,tmp);
            }
        }
    }

    for (uint_t l = 0;l<settings.levels;l++){
        delete domains[l];
    }

    cout << "Maximal error: " << err << endl;

    return EXIT_SUCCESS;
}