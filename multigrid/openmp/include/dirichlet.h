#ifndef DIRICHLET_BOUNDARY
#define DIRICHLET_BOUNDARY

#include "boundary.h"

template <class T>
class Dirichlet :
    public Boundary<T> {
        public:
        Dirichlet(int_t device,Location_t location,uint_t i, uint_t j, uint_t k);

        virtual ~Dirichlet();

        void init(funptr ufun,funptr dudxfun,funptr dudyfun,Settings & settings);

};

template <class T>
Dirichlet<T>::Dirichlet(int_t device,Location_t location,uint_t i, uint_t j, uint_t k) 
    : Boundary<T>(device,location,i,j,k){

}

template <class T>
Dirichlet<T>::~Dirichlet(){
    // Does nothing
}

template <class T>
void Dirichlet<T>::init(funptr ufun,funptr dudxfun,funptr dudyfun,Settings & settings){
    this->init_by_fun(ufun,settings);
}

#endif
