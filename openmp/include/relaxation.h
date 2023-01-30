#ifndef POISSON_RELAXATION
#define POISSON_RELAXATION

#include "domain.h"

namespace Poisson{
    template <class T>
    class Relaxation {
            public:
            virtual void relax(Domain<T>& domain,T omega);
            virtual constexpr T default_omega()=0;
            virtual constexpr bool requires_duplicate_solution()=0;
    };
}

#endif