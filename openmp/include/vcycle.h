#ifndef POISSON_VCYCLE
#define POISSON_VCYCLE

#include "residual.h"
#include "jacobi.h"
#include "gaussseidel.h"
#include <iostream>

using std::cout;
using std::endl;

namespace Poisson{
    template <class T>
    void Vcycle(
        Domain<T> *** domains,
        Restriction<T> & restriction,
        Prolongation<T> & prolongation,
        Relaxation<T> & relaxation,
        T omega,
        uint_t level,
        uint_t levels,
        uint_t num_devices,
        int_t nsmooth = 10){
        
        // Pre smooting
        
        relaxation.relax(domains,omega,nsmooth,level,num_devices);

        if (level >= levels-1){
            return;
        }

        // Compute and restrict the defect
        #pragma omp parallel for num_threads(num_devices)
        for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
            residual<T>(*domains[gpuid][level]);
            
            restriction.restrict_to(*(domains[gpuid][level]->r),*(domains[gpuid][level+1]->f));
            
            domains[gpuid][level+1]->u->init_zero();
            
            domains[gpuid][level+1]->north->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->north),(*domains[gpuid][level]).settings,restriction);
            domains[gpuid][level+1]->south->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->south),(*domains[gpuid][level]).settings,restriction);
            domains[gpuid][level+1]->east->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->east),(*domains[gpuid][level]).settings,restriction);
            domains[gpuid][level+1]->west->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->west),(*domains[gpuid][level]).settings,restriction);
            domains[gpuid][level+1]->top->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->top),(*domains[gpuid][level]).settings,restriction);
            domains[gpuid][level+1]->bottom->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->bottom),(*domains[gpuid][level]).settings,restriction);
        }
        // Recursion
        Vcycle<T>(domains,restriction,prolongation,relaxation,omega,level+1,levels,num_devices,nsmooth);

        // Interpolate error
        #pragma omp parallel for num_threads(num_devices)
        for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
            prolongation.prolong(*(domains[gpuid][level+1]->u),*(domains[gpuid][level]->r));
            domains[gpuid][level]->u->add(*(domains[gpuid][level]->r));
        }

        // Post smooting
        relaxation.relax(domains,omega,nsmooth,level,num_devices);
    }
}

#endif
