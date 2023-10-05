#ifndef DOMAIN_HALO
#define DOMAIN_HALO

#include "definitions.h"
namespace Poisson{
    class Halo {
        public:
            const uint_t east;
            const uint_t west;
            const uint_t north;
            const uint_t south;
            const uint_t top;
            const uint_t bottom;

            Halo() : 
            east(0), west(0),  north(0), south(0), top(0), bottom(0) {

            }
        
            Halo(uint_t _east, uint_t _west, uint_t _north, uint_t _south, uint_t _top,uint_t _bottom) : 
            east(_east), west(_west),  north(_north), south(_south), top(_top), bottom(_bottom) {

            }
            // To ensure that the class is trivially copyable, the copy constructor has been commented out 
            // Halo(Halo & _halo) : 
            //     east(_halo.east), west(_halo.west), north(_halo.north), south(_halo.south), top(_halo.top), bottom(_halo.bottom) {

            // }

            // To ensure that the class is trivially copyable, the destructor has been commented out
            //~Halo(){

            //}
    };
}

#endif
