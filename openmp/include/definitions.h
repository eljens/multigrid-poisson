#ifndef POISSON_DEFINITIONS
#define POISSON_DEFINITIONS

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include <cstdint>
namespace Poisson{
    #ifdef SCHED
    #define SCHEDULE schedule(SCHED)
    #else
    #define SCHEDULE /**/
    #endif

    #ifdef DIST
    #define DIST_SCHEDULE dist_schedule(DIST)
    #else
    #define DIST_SCHEDULE /**/
    #endif

    using double_t = double;
    using float_t = float;
    using uint_t = std::size_t;
    using int_t = long int;
}
#endif
