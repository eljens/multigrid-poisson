#ifndef POISSON_BOUNDARY
#define POISSON_BOUNDARY

#include "array.h"

typedef enum {NORTH,SOUTH,EAST,WEST,TOP,BOTTOM} Location_t;

template <class T>
class Boundary {
    protected:
    Array<T> * at;
    Location_t location;
    public:
    Boundary(Location_t location,uint_t dim1,uint_t dim2){
        this->location = location;
        at = new Array<T>({dim1,dim2});
    }
    ~Boundary(){
        delete at;
    }
    //virtual void write() = 0;
};

#endif