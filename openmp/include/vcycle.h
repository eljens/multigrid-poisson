#ifndef POISSON_VCYCLE
#define POISSON_VCYCLE

#include "residual.h"
#include "jacobi.h"
#include <iostream>

#define nsmooth 5

using std::cout;
using std::endl;

template <class T>
void Vcycle(
    Domain<T> * domains[],
    Restriction<T> & restriction,
    Prolongation<T> & prolongation,
    T omega,
    uint_t level,
    uint_t levels){
    
    // Pre smooting
    const uint_t _nsmooth = (level == levels-1) ? nsmooth : nsmooth;
    for(int_t i=0;i<_nsmooth;i++){
        jacobi<double_t>(*domains[level],omega);
    }

    if (level >= levels-1){
        return;
    }

    residual<T>(*domains[level]);

    restriction.restrict(*(domains[level]->r),*(domains[level+1]->f));

    domains[level+1]->u->init_zero();
    domains[level+1]->uprev->init_zero();

    // Restricting boundaries
    domains[level+1]->north->restrict(*(domains[level]->u),*(domains[level]->north),(*domains[level]).settings,restriction);
    domains[level+1]->south->restrict(*(domains[level]->u),*(domains[level]->south),(*domains[level]).settings,restriction);
    domains[level+1]->east->restrict(*(domains[level]->u),*(domains[level]->east),(*domains[level]).settings,restriction);
    domains[level+1]->west->restrict(*(domains[level]->u),*(domains[level]->west),(*domains[level]).settings,restriction);
    domains[level+1]->top->restrict(*(domains[level]->u),*(domains[level]->top),(*domains[level]).settings,restriction);
    domains[level+1]->bottom->restrict(*(domains[level]->u),*(domains[level]->bottom),(*domains[level]).settings,restriction);

    // Recursion
    Vcycle<T>(domains,restriction,prolongation,omega,level+1,levels);

    prolongation.prolong(*(domains[level+1]->u),*(domains[level]->r));

    // Interpolate error
    domains[level]->uprev->add(*(domains[level]->r));

    // Post smooting
    for(int_t i=0;i<nsmooth;i++){
        jacobi<double_t>(*domains[level],omega);
    }
}

#endif