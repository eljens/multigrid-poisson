#ifndef POISSON_FULL_WEIGHT
#define POISSON_INJECTION

#include "restriction.h"

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

namespace Poisson{
    template<class T>
    class FullWeighting :
        public Restriction<T> {
            public:
            FullWeighting();

            ~FullWeighting();

            void restrict_to(DeviceArray<T> & uin,DeviceArray<T> & uout);
    };

    template<class T>
    FullWeighting<T>::FullWeighting() : Restriction<T>() {
        // Does nothing
    }

    template<class T>
    FullWeighting<T>::~FullWeighting(){
        // Does nothing
    }

    template<class T>
    void FullWeighting<T>::restrict_to(DeviceArray<T> & uin,DeviceArray<T> & uout){
        bool cond0 = ((uint_t) (uin.shape[0]-1)/2)+1 == uout.shape[0];
        bool cond1 = ((uint_t) (uin.shape[1]-1)/2)+1 == uout.shape[1];
        bool cond2 = ((uint_t) (uin.shape[2]-1)/2)+1 == uout.shape[2];
        if ((!cond0) || (!cond1) || (!cond2)){
            cerr << "Restriction from domain of size (" << uin.shape[0] << "," << uin.shape[1] << "," << uin.shape[2] << ")";
            cerr << " to (" << uout.shape[0] << "," << uout.shape[1] << "," << uout.shape[2] << ") is incompatible with full weighting" << endl;
        }

        T * indev = uin.devptr;
        T * outdev = uout.devptr;

        // Dimensions of output
        const uint_t (&uoutshape)[3] = uout.shape;
        const Halo & uouthalo = uout.halo;
        const uint_t (&uoutstride)[3] = uout.stride;

        // Dimensions of input
        const uint_t (&uinshape)[3] = uin.shape;
        const Halo & uinhalo = uin.halo;
        const uint_t (&uinstride)[3] = uin.stride;

        #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE) dist_schedule(static,DIST_SIZE)
            for (int_t i = 0;i<uoutshape[0];i++){
                for (int_t j = 0;j<uoutshape[1];j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 0;k_block<uoutshape[2];k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uoutshape[2]);k++){
#else
                    for (int_t k = 0;k<uoutshape[2];k++){
#endif
                            int_t ii = MIN(uinshape[0]-1,2*i);
                            int_t jj = MIN(uinshape[1]-1,2*j);
                            int_t kk = MIN(uinshape[2]-1,2*k);
                            int_t ii_m1 = MAX(0,ii-1);
                            int_t ii_p1 = MIN(uinshape[0]-1,ii+1);
                            int_t jj_m1 = MAX(0,jj-1);
                            int_t jj_p1 = MIN(uinshape[1]-1,jj+1);
                            int_t kk_m1 = MAX(0,kk-1);
                            int_t kk_p1 = MIN(uinshape[2]-1,kk+1);

                            // -1 offset in z
                            T tmp = 1.0*indev[idx(ii_m1,jj_m1,kk_m1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_m1,jj,kk_m1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_m1,jj_p1,kk_m1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii,jj_m1,kk_m1,uinhalo,uinstride)] +
                                    4.0*indev[idx(ii,jj,kk_m1,uinhalo,uinstride)]+
                                    2.0*indev[idx(ii,jj_p1,kk_m1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_p1,jj_m1,kk_m1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_p1,jj,kk_m1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_p1,jj_p1,kk_m1,uinhalo,uinstride)];

                            // 0 offset in z
                            tmp +=  2.0*indev[idx(ii_m1,jj_m1,kk,uinhalo,uinstride)] +
                                    4.0*indev[idx(ii_m1,jj,kk,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_m1,jj_p1,kk,uinhalo,uinstride)] +
                                    4.0*indev[idx(ii,jj_m1,kk,uinhalo,uinstride)] +
                                    8.0*indev[idx(ii,jj,kk,uinhalo,uinstride)]+
                                    4.0*indev[idx(ii,jj_p1,kk,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_p1,jj_m1,kk,uinhalo,uinstride)] +
                                    4.0*indev[idx(ii_p1,jj,kk,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_p1,jj_p1,kk,uinhalo,uinstride)];
                            
                            // +1 offset in z
                            tmp +=  1.0*indev[idx(ii_m1,jj_m1,kk_p1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_m1,jj,kk_p1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_m1,jj_p1,kk_p1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii,jj_m1,kk_p1,uinhalo,uinstride)] +
                                    4.0*indev[idx(ii,jj,kk_p1,uinhalo,uinstride)]+
                                    2.0*indev[idx(ii,jj_p1,kk_p1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_p1,jj_m1,kk_p1,uinhalo,uinstride)] +
                                    2.0*indev[idx(ii_p1,jj,kk_p1,uinhalo,uinstride)] +
                                    1.0*indev[idx(ii_p1,jj_p1,kk_p1,uinhalo,uinstride)];

                            outdev[idx(i,j,k,uouthalo,uoutstride)] = tmp/64.0;
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
