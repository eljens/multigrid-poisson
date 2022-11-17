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

        void write_to(DeviceArray<T> & uarr, Settings & settings);

        void update(DeviceArray<T> & uarr, Settings & settings);
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
            cerr << "Neumann conditions are not yet supported on TOP and BOTTOM" << endl; 
            break;
    }
}

template<class T>
void Neumann<T>::write_to(DeviceArray<T> & uarr, Settings & settings){
    int_t offx = uarr.halo.west;
    int_t offy = uarr.halo.south;
    int_t offz = uarr.halo.bottom;

    int_t ii = 0;
    int_t jj = 0;
    int_t kk = 0;

    T sign = 1.0;

    switch (this->location){
        case EAST:
            offx = uarr.shape[0]-1+uarr.halo.east+uarr.halo.west;
            ii = -2;
            break;
        case WEST:
            offx = 0;
            ii = 2;
            sign = -1.0;
            break;
        case NORTH:
            offy = uarr.shape[1]-1+uarr.halo.north+uarr.halo.south;
            jj = -2;
            break;
        case SOUTH:
            offy=0;
            sign = -1.0;
            jj = 2;
            break;
        case TOP:
            offz = uarr.shape[2]-1+uarr.halo.top+uarr.halo.bottom;
            kk = -2;
            break;
        case BOTTOM:
            offz = 0;
            sign = -1.0;
            kk = 2;
            break;
    }
    T * udev = uarr.devptr;
    T * gdev = this->devptr;
    const T two_h = 2.0*settings.h;
    #pragma omp target device(this->device) is_device_ptr(udev,gdev) firstprivate(sign,two_h)
    #pragma omp teams distribute parallel for collapse(3) schedule(static,1)
    for(int_t i = 0;i<this->shape[0];i++){
        for(int_t j = 0;j<this->shape[1];j++){
            for(int_t k = 0;k<this->shape[2];k++){
                udev[uarr.idx_halo(i+offx,j+offy,k+offz)] =
                    sign*two_h*gdev[this->idx(i,j,k)]+udev[uarr.idx_halo(i+offx+ii,j+offy+jj,k+offz+kk)];
            }
        }
    }
};

template<class T>
void Neumann<T>::update(DeviceArray<T> & uarr, Settings & settings){
    this->write_to(uarr,settings);
}

#endif