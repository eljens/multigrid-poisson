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
        for(int_t i=0;i<nsmooth;i++){
            #pragma omp parallel
            #pragma omp single nowait
            {
                for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
                    #pragma omp task
                    relaxation.relax(*domains[gpuid][level],omega);
                }
            }
        }

        if (level >= levels-1){
            return;
        }

        // Compute and restrict the defect
        #pragma omp parallel
        #pragma omp single nowait
        {
            for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
                #pragma omp task default(none) shared(domains,level) firstprivate(gpuid) depend(inout:domains[gpuid][level]->r)
                {
                    residual<T>(*domains[gpuid][level]);
                }
                #pragma omp task default(none) shared(restriction,domains,level) firstprivate(gpuid) depend(in:domains[gpuid][level]->r) depend(out:domains[gpuid][level+1]->f)
                {
                    restriction.restrict_to(*(domains[gpuid][level]->r),*(domains[gpuid][level+1]->f));
                }

                //domains[level+1]->u->init_zero();
                #pragma omp task
                domains[gpuid][level+1]->u->init_zero();

                // Restricting boundaries
                #pragma omp task
                domains[gpuid][level+1]->north->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->north),(*domains[gpuid][level]).settings,restriction);
                #pragma omp task
                domains[gpuid][level+1]->south->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->south),(*domains[gpuid][level]).settings,restriction);
                #pragma omp task
                domains[gpuid][level+1]->east->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->east),(*domains[gpuid][level]).settings,restriction);
                #pragma omp task
                domains[gpuid][level+1]->west->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->west),(*domains[gpuid][level]).settings,restriction);
                #pragma omp task
                domains[gpuid][level+1]->top->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->top),(*domains[gpuid][level]).settings,restriction);
                #pragma omp task
                domains[gpuid][level+1]->bottom->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->bottom),(*domains[gpuid][level]).settings,restriction);
            }
        }
        // Recursion
        Vcycle<T>(domains,restriction,prolongation,relaxation,omega,level+1,levels,num_devices,nsmooth);

        // Interpolate error
        #pragma omp parallel
        #pragma omp single nowait
        {
            for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
                #pragma omp task default(none) shared(prolongation,domains,level) firstprivate(gpuid) depend(in:domains[gpuid][level+1]->u) depend(out:domains[gpuid][level]->r)
                {
                    prolongation.prolong(*(domains[gpuid][level+1]->u),*(domains[gpuid][level]->r));
                }
                #pragma omp task default(none) shared(domains,level) firstprivate(gpuid) depend(in:domains[gpuid][level]->r) depend(out:domains[gpuid][level]->u)
                {
                    domains[gpuid][level]->u->add(*(domains[gpuid][level]->r));
                }
            }
        }

        // Post smooting
        for(int_t i=0;i<nsmooth;i++){
            #pragma omp parallel
            #pragma omp single nowait
            {
                for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
                    #pragma omp task
                    relaxation.relax(*domains[gpuid][level],omega);
                }
            }
        }
    }
}

#endif
