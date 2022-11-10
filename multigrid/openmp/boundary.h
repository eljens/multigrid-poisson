#ifndef POISSON_BOUNDARY
#define POISSON_BOUNDARY

#include "array.h"

typedef enum {NORTH,SOUTH,EAST,WEST,TOP,BOTTOM} Location_t;

template <class T>
class Boundary :
    public Array<T>{
    protected:
        Location_t location;
        public:
        Boundary(Location_t location,initializer_list<uint_t> args,Array<T> domain) : Array<T>(args){
            this->location = location;
        }
        ~Boundary(){
            delete at;
        }
        //virtual void write() = 0;
};

#endif