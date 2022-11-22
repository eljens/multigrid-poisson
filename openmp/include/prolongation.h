#ifndef POISSON_PROLONGATION
#define POISSON_PROLONGATION

#include "devicearray.h"

template <class T>
class Prolongation {
        public:
        Prolongation();

        virtual ~Prolongation();

        virtual void prolong(DeviceArray<T> & uin,DeviceArray<T> & uout) = 0;
};

template<class T>
Prolongation<T>::Prolongation() {
    // Does nothing
}

template<class T>
Prolongation<T>::~Prolongation() {
    // Does nothing
}

#endif