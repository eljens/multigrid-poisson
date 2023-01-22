#ifndef POISSON_TRILINEAR_INTERPOLATION
#define POISSON_TRILINEAR_INTERPOLATION

#include "prolongation.h"

namespace Poisson{
    template<class T>
    class TrilinearInterpolation :
        public Prolongation<T> {
            public:
            TrilinearInterpolation();

            ~TrilinearInterpolation();

            void prolong(DeviceArray<T> & uin,DeviceArray<T> & uout);
    };

    template<class T>
    TrilinearInterpolation<T>::TrilinearInterpolation() : Prolongation<T>() {
        // Does nothing
    }

    template<class T>
    TrilinearInterpolation<T>::~TrilinearInterpolation(){
        // Does nothing
    }

    template<class T>
    void TrilinearInterpolation<T>::prolong(DeviceArray<T> & uin,DeviceArray<T> & uout){
        bool cond0 = ((uint_t) (uout.shape[0]-1)/2)+1 == uin.shape[0];
        bool cond1 = ((uint_t) (uout.shape[1]-1)/2)+1 == uin.shape[1];
        bool cond2 = ((uint_t) (uout.shape[2]-1)/2)+1 == uin.shape[2];
        if ((!cond0) || (!cond1) || (!cond2)){
            cerr << "Prolongation from domain of size (" << uin.shape[0] << "," << uin.shape[1] << "," << uin.shape[2] << ")";
            cerr << " to (" << uout.shape[0] << "," << uout.shape[1] << "," << uout.shape[2] << ") is incompatible with trilinear interpolation" << endl;
        }

        T * indev = uin.devptr;
        T * outdev = uout.devptr;

        #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
        {
            //uf(1:2:end,1:2:end,1:2:end) = uc;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uin.shape[0];i++){
                for (int_t j = 0;j<uin.shape[1];j++){
                    for (int_t k = 0;k<uin.shape[2];k++){
                        outdev[uout.idx(2*i,2*j,2*k)] = indev[uin.idx(i,j,k)];
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(2:2:end,:,:) = (uf(3:2:end,:,:) + uf(1:2:end-2,:,:))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 1;i<uout.shape[0]-1;i+=2){
                for (int_t j = 0;j<uout.shape[1];j+=2){
                    for (int_t k = 0;k<uout.shape[2];k+=2){
                        outdev[uout.idx(i,j,k)] = 0.5*(outdev[uout.idx(i-1,j,k)]+outdev[uout.idx(i+1,j,k)]);
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(:,2:2:end,:) = (uf(:,3:2:end,:) + uf(:,1:2:end-2,:))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uout.shape[0];i++){
                for (int_t j = 1;j<uout.shape[1]-1;j+=2){
                    for (int_t k = 0;k<uout.shape[2];k+=2){
                        outdev[uout.idx(i,j,k)] = 0.5*(outdev[uout.idx(i,j-1,k)]+outdev[uout.idx(i,j+1,k)]);
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(:,:,2:2:end) = (uf(:,:,3:2:end) + uf(:,:,1:2:end-2))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uout.shape[0];i++){
                for (int_t j = 0;j<uout.shape[1];j++){
                    for (int_t k = 1;k<uout.shape[2]-1;k+=2){
                        outdev[uout.idx(i,j,k)] = 0.5*(outdev[uout.idx(i,j,k-1)]+outdev[uout.idx(i,j,k+1)]);
                    }
                }
            }
        }
    }
}
#endif