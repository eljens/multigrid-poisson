#ifndef POISSON_FCYCLE
#define POISSON_FCYCLE

#include "residual.h"
#include "relaxation.h"
#include <iostream>

using std::cout;
using std::endl;

namespace Poisson{
    template <class T>
    void Fcycle(
        Domain<T> *** domains,
        Restriction<T> & restriction,
        Prolongation<T> & prolongation,
        Relaxation<T> & relaxation,
        T omega,
        uint_t level,
        uint_t levels,
        uint_t num_devices,
        int_t nsmooth = 10)
    {    

        // Pre smooting
        relaxation.relax(domains,omega,nsmooth,level,num_devices);

        // Compute and restrict the defect
        #pragma omp parallel
        #pragma omp single nowait
        #pragma omp taskgroup
        for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
            #pragma omp task default(none) shared(domains,level) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level]->r->at[0])
            residual<T>(*domains[gpuid][level]);
            
            #pragma omp task default(none) shared(restriction,domains,level) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->r->at[0])\
            depend(out:domains[gpuid][level+1]->f->at[0])
            restriction.restrict_to(*(domains[gpuid][level]->r),*(domains[gpuid][level+1]->f));
            
            #pragma omp task default(none) shared(domains,level) firstprivate(gpuid)\
            depend(out:domains[gpuid][level+1]->u->at[0])
            domains[gpuid][level+1]->u->init_zero();
            
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->north->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->north->arr.at[0])
            domains[gpuid][level+1]->north->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->north),(*domains[gpuid][level]).settings,restriction);
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->south->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->south->arr.at[0])
            domains[gpuid][level+1]->south->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->south),(*domains[gpuid][level]).settings,restriction);
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->east->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->east->arr.at[0])
            domains[gpuid][level+1]->east->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->east),(*domains[gpuid][level]).settings,restriction);
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->west->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->west->arr.at[0])
            domains[gpuid][level+1]->west->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->west),(*domains[gpuid][level]).settings,restriction);
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->top->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->top->arr.at[0])
            domains[gpuid][level+1]->top->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->top),(*domains[gpuid][level]).settings,restriction);
            #pragma omp task default(none) shared(domains,level,restriction) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->bottom->arr.at[0],domains[gpuid][level]->u->at[0])\
            depend(out:domains[gpuid][level+1]->bottom->arr.at[0])
            domains[gpuid][level+1]->bottom->restrict_to(*(domains[gpuid][level]->u),*(domains[gpuid][level]->bottom),(*domains[gpuid][level]).settings,restriction);
        }
        if (level < levels-2){
            Fcycle(domains,restriction,prolongation,relaxation,omega,level+1,levels,num_devices,nsmooth);
        }

        // Recursion
        Vcycle<T>(domains,restriction,prolongation,relaxation,omega,level+1,levels,num_devices,nsmooth);

                // Interpolate error
        #pragma omp parallel
        #pragma omp single nowait
        #pragma omp taskgroup
        for (int_t gpuid = 0; gpuid < num_devices;gpuid++){
            #pragma omp task default(none) shared(prolongation,domains,level) firstprivate(gpuid)\
            depend(in:domains[gpuid][level+1]->u->at[0])\
            depend(out:domains[gpuid][level]->r->at[0])
            prolongation.prolong(*(domains[gpuid][level+1]->u),*(domains[gpuid][level]->r));
            
            #pragma omp task default(none) shared(domains,level) firstprivate(gpuid)\
            depend(in:domains[gpuid][level]->r->at[0])\
            depend(out:domains[gpuid][level]->u->at[0])
            domains[gpuid][level]->u->add(*(domains[gpuid][level]->r));
        }
        if (level == 0){
            Vcycle<T>(domains,restriction,prolongation,relaxation,omega,0,levels,num_devices,nsmooth);
        }
    }
}

#endif