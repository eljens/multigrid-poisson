#ifndef DOMAIN_HALO
#define DOMAIN_HALO

#include "definitions.h"

class Halo {
    public:
        const uint_t east;
        const uint_t west;
        const uint_t north;
        const uint_t south;
        const uint_t top;
        const uint_t bottom;
    
        Halo(uint_t _east, uint_t _west, uint_t _north, uint_t _south, uint_t _top,uint_t _bottom) : 
           east(_east), west(_west),  north(_north), south(_south), top(_top), bottom(_bottom) {

        }

        Halo(Halo & _halo) : 
            east(_halo.east), west(_halo.west), north(_halo.north), south(_halo.south), top(_halo.top), bottom(_halo.bottom) {

        }

        ~Halo(){

        }
};

#endif