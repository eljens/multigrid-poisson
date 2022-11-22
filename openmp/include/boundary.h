#ifndef POISSON_BOUNDARY
#define POISSON_BOUNDARY

#include "devicearray.h"
#include "definitions.h"
#include "problem_definition.h"
#include "settings.h"
#include "halo.h"
#include "restriction.h"

typedef enum {NORTH,SOUTH,EAST,WEST,TOP,BOTTOM} Location_t;

template <class T>
class Boundary :
    public DeviceArray<T> {
    protected:
        Location_t location;
        void init_by_fun(funptr fun,Settings & settings);

        public:
        Boundary(int_t device, Location_t location,uint_t i, uint_t j, uint_t k);

        virtual ~Boundary();
        
        virtual void init(funptr ufun,funptr dudxfun,funptr dudyfun,Settings & settings);

        virtual void write_to(DeviceArray<T> & uarr, Settings & settings);

        virtual void update(DeviceArray<T> & uarr, Settings & settings);

        virtual void restrict(DeviceArray<T> & u, Boundary<T> & boundary,
                                Settings & settings, Restriction<T> & restriction);
};

template<class T>
Boundary<T>::Boundary(int_t device, Location_t location,uint_t i, uint_t j, uint_t k) :
    DeviceArray<T>(device,i,j,k) {
    this->location = location;
}

template<class T>
Boundary<T>::~Boundary() {
    // Does nothing
}

template<class T>
void Boundary<T>::init_by_fun(funptr fun,Settings & settings){
    const double_t x0 = settings.origin[0];
    const double_t y0 = settings.origin[1];
    const double_t z0 = settings.origin[2];
    const double_t h = settings.h;
    double_t offx = 0;
    double_t offy = 0;
    double_t offz = 0;
    switch (this->location){
        case TOP:
            offz = settings.dims[2]-1;
            break;
        case NORTH:
            offy = settings.dims[1]-1;
            break;
        case EAST:
            offx = settings.dims[0]-1;
            break;
        default:
            break;
    }
    for(int_t i = 0;i<this->shape[0];i++){
        for(int_t j = 0;j<this->shape[1];j++){
            for(int_t k = 0;k<this->shape[2];k++){
                this->at[this->idx(i,j,k)] = fun(x0+(((double_t)i)+offx)*h,y0+(((double_t)j)+offy)*h,z0+(((double_t)k)+offz)*h);
            }
        }
    }
}

#endif