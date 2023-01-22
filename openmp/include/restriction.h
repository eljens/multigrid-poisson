#ifndef POISSON_RESTRICTION
#define POISSON_RESTRICTION

#include "devicearray.h"
namespace Poisson{
    template <class T>
    class Restriction {
            public:
            Restriction();

            virtual ~Restriction();

            virtual void restrict(DeviceArray<T> & uin,DeviceArray<T> & uout) = 0;
    };

    template<class T>
    Restriction<T>::Restriction() {
        // Does nothing
    }

    template<class T>
    Restriction<T>::~Restriction() {
        // Does nothing
    }
}
#endif