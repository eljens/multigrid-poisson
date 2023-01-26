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
            cerr << " to (" << uout.shape[0] << "," << uout.shape[1] << "," << uout.shape[2] << ") is incompatible with injection" << endl;
        }

        T * indev = uin.devptr;
        T * outdev = uout.devptr;

        #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uout.shape[0];i++){
                for (int_t j = 0;j<uout.shape[1];j++){
                    for (int_t k = 0;k<uout.shape[2];k++){
                        int_t ii = 2*i;
                        int_t jj = 2*j;
                        int_t kk = 2*k;
                        int_t ii_m1 = MAX(0,ii-1);
                        int_t ii_p1 = MIN(uin.shape[0]-1,ii+1);
                        int_t jj_m1 = MAX(0,jj-1);
                        int_t jj_p1 = MIN(uin.shape[1]-1,jj+1);
                        int_t kk_m1 = MAX(0,kk-1);
                        int_t kk_p1 = MIN(uin.shape[2]-1,kk+1);

                        // -1 offset in z
                        T tmp = 1.0*indev[uin.idx(ii_m1,jj_m1,kk_m1)] +
                                2.0*indev[uin.idx(ii_m1,jj,kk_m1)] +
                                1.0*indev[uin.idx(ii_m1,jj_p1,kk_m1)] +
                                2.0*indev[uin.idx(ii,jj_m1,kk_m1)] +
                                4.0*indev[uin.idx(ii,jj,kk_m1)]+
                                2.0*indev[uin.idx(ii,jj_p1,kk_m1)] +
                                1.0*indev[uin.idx(ii_p1,jj_m1,kk_m1)] +
                                2.0*indev[uin.idx(ii_p1,jj,kk_m1)] +
                                1.0*indev[uin.idx(ii_p1,jj_p1,kk_m1)];

                        // 0 offset in z
                        tmp +=  2.0*indev[uin.idx(ii_m1,jj_m1,kk)] +
                                4.0*indev[uin.idx(ii_m1,jj,kk)] +
                                2.0*indev[uin.idx(ii_m1,jj_p1,kk)] +
                                4.0*indev[uin.idx(ii,jj_m1,kk)] +
                                8.0*indev[uin.idx(ii,jj,kk)]+
                                4.0*indev[uin.idx(ii,jj_p1,kk)] +
                                2.0*indev[uin.idx(ii_p1,jj_m1,kk)] +
                                4.0*indev[uin.idx(ii_p1,jj,kk)] +
                                2.0*indev[uin.idx(ii_p1,jj_p1,kk)];
                        
                        // +1 offset in z
                        tmp +=  1.0*indev[uin.idx(ii_m1,jj_m1,kk_p1)] +
                                2.0*indev[uin.idx(ii_m1,jj,kk_p1)] +
                                1.0*indev[uin.idx(ii_m1,jj_p1,kk_p1)] +
                                2.0*indev[uin.idx(ii,jj_m1,kk_p1)] +
                                4.0*indev[uin.idx(ii,jj,kk_p1)]+
                                2.0*indev[uin.idx(ii,jj_p1,kk_p1)] +
                                1.0*indev[uin.idx(ii_p1,jj_m1,kk_p1)] +
                                2.0*indev[uin.idx(ii_p1,jj,kk_p1)] +
                                1.0*indev[uin.idx(ii_p1,jj_p1,kk_p1)];

                        outdev[uout.idx(i,j,k)] = tmp/64.0;
                    }
                }
            }
        }
    }
}
#endif
