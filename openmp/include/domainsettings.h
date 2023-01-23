#ifndef DOMAIN_SETTINGS
#define DOMAIN_SETTINGS

#include <iostream>
#include <string>
#include "definitions.h"
#include "settings.h"

namespace Poisson{
    using std::cout;
    using std::endl;
    using std::ostream;

    class DomainSettings :
        public Settings {
            public:
                DomainSettings();

                DomainSettings(Settings & settings, const uint_t l);

                ~DomainSettings();
    };

    ostream& operator<<(ostream& os, const DomainSettings& settings);
}
#endif
