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
        Domain<T> * domains[],
        Restriction<T> & restriction,
        Prolongation<T> & prolongation,
        Relaxation<T> & relaxation,
        T omega,
        uint_t level,
        uint_t levels,
        int_t nsmooth = 10){
        
        // Pre smooting
        for(int_t i=0;i<nsmooth;i++){
            relaxation.relax(*domains[level],omega);
        }

        if (level >= levels-1){
            return;
        }

        // Compute and restrict the defect
        residual<T>(*domains[level]);
        restriction.restrict_to(*(domains[level]->r),*(domains[level+1]->f));

        //domains[level+1]->u->init_zero();
        domains[level+1]->u->init_zero();

        // Restricting boundaries
        domains[level+1]->north->restrict_to(*(domains[level]->u),*(domains[level]->north),(*domains[level]).settings,restriction);
        domains[level+1]->south->restrict_to(*(domains[level]->u),*(domains[level]->south),(*domains[level]).settings,restriction);
        domains[level+1]->east->restrict_to(*(domains[level]->u),*(domains[level]->east),(*domains[level]).settings,restriction);
        domains[level+1]->west->restrict_to(*(domains[level]->u),*(domains[level]->west),(*domains[level]).settings,restriction);
        domains[level+1]->top->restrict_to(*(domains[level]->u),*(domains[level]->top),(*domains[level]).settings,restriction);
        domains[level+1]->bottom->restrict_to(*(domains[level]->u),*(domains[level]->bottom),(*domains[level]).settings,restriction);

        // Recursion
        Vcycle<T>(domains,restriction,prolongation,relaxation,omega,level+1,levels);

        // Interpolate error
        prolongation.prolong(*(domains[level+1]->u),*(domains[level]->r));
        domains[level]->u->add(*(domains[level]->r));

        // Post smooting
        for(int_t i=0;i<nsmooth;i++){
            relaxation.relax(*domains[level],omega);
        }
    }
}

#endif
