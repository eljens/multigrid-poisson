#ifndef POISSON_DEFINITIONS
#define POISSON_DEFINITIONS

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include <cstdint>
namespace Poisson{
    #ifndef CHUNK_SIZE
    #define CHUNK_SIZE 1
    #endif

    #ifdef DIST_SIZE
    #define DIST_SCHEDULE dist_schedule(static,DIST_SIZE)
    #else
    #define DIST_SCHEDULE /**/
    #endif

    using double_t = double;
    using float_t = float;
    using uint_t = std::size_t;
    using int_t = long int;
}
#endif
