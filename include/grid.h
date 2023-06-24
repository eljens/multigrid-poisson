#ifndef POISSON_GRID
#define POISSON_GRID

#include "settings.h"
#include "domainsettings.h"

namespace Poisson{
    class Grid {
        public:
            const uint_t levels;
            DomainSettings * domainsettings;

            Grid(const Settings & settings,const uint_t _levels);

            ~Grid();
    };
}

#endif
