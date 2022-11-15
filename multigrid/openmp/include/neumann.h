#ifndef NEUMANN_BOUNDARY
#define NEUMANN_BOUNDARY

#include "boundary.h"

template <class T>
class Neumann :
    public Boundary<T> {
        public:
        Neumann(int_t device,Location_t location,uint_t i, uint_t j, uint_t k);

        virtual ~Neumann();

        void init(funptr ufun,funptr dudxfun,funptr dudyfun,Settings & settings);
};

template<class T>
Neumann<T>::Neumann(int_t device,Location_t location,uint_t i, uint_t j, uint_t k) 
    : Boundary<T>(device,location,i,j,k){

};

template<class T>
Neumann<T>::~Neumann() {
    // Does nothing
}

template<class T>
void Neumann<T>::init(funptr ufun,funptr dudxfun,funptr dudyfun,Settings & settings){
    switch (this->location){
        case EAST:
        case WEST:
            this->init_by_fun(dudxfun,settings);
            break;
        case NORTH:
        case SOUTH:
            this->init_by_fun(dudyfun,settings);
            break;
        default:
            break;
    }
}

#endif