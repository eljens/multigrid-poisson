#ifndef POISSON_INJECTION
#define POISSON_INJECTION

#include "restriction.h"

namespace Poisson{
    template<class T>
    class Injection :
        public Restriction<T> {
            public:
            Injection();

            ~Injection();

            void restrict_to(DeviceArray<T> & uin,DeviceArray<T> & uout);
    };

    template<class T>
    Injection<T>::Injection() : Restriction<T>() {
        // Does nothing
    }

    template<class T>
    Injection<T>::~Injection(){
        // Does nothing
    }

    template<class T>
    void Injection<T>::restrict_to(DeviceArray<T> & uin,DeviceArray<T> & uout){
        bool cond0 = ((uint_t) (uin.shape[0]-1)/2)+1 == uout.shape[0];
        bool cond1 = ((uint_t) (uin.shape[1]-1)/2)+1 == uout.shape[1];
        bool cond2 = ((uint_t) (uin.shape[2]-1)/2)+1 == uout.shape[2];
        if ((!cond0) || (!cond1) || (!cond2)){
            cerr << "Restriction from domain of size (" << uin.shape[0] << "," << uin.shape[1] << "," << uin.shape[2] << ")";
            cerr << " to (" << uout.shape[0] << "," << uout.shape[1] << "," << uout.shape[2] << ") is incompatible with injection" << endl;
        }

        T * indev = uin.devptr;
        T * outdev = uout.devptr;
        // Shape of output array
        const uint_t (&uoutshape)[3] = uout.shape;
        const Halo & uouthalo = uout.halo;
        const uint_t (&uoutstride)[3] = uout.stride;
        // Shape of input array
        const Halo & uinhalo = uin.halo;
        const uint_t (&uinstride)[3] = uin.stride;
        #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE) DIST_SCHEDULE
            for (int_t i = 0;i<uoutshape[0];i++){
                for (int_t j = 0;j<uoutshape[1];j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 0;k_block<uoutshape[2];k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uoutshape[2]);k++){
#else
                    for (int_t k = 0;k<uoutshape[2];k++){
#endif
                            outdev[idx(i,j,k,uouthalo,uoutstride)] = indev[idx(2*i,2*j,2*k,uinhalo,uinstride)];
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
    }
}
#endif
