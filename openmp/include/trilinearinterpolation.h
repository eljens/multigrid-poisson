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

        // Shape, strid and halo size for arrays
        const Halo & uinhalo = uin.halo;
        const uint_t (&uinstride)[3] = uin.stride;
        const uint_t (&uinshape)[3] = uin.shape;

        const Halo & uouthalo = uout.halo;
        const uint_t (&uoutstride)[3] = uout.stride;
        const uint_t (&uoutshape)[3] = uout.shape;

        #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
        {
            //uf(1:2:end,1:2:end,1:2:end) = uc;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uinshape[0];i++){
                for (int_t j = 0;j<uinshape[1];j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 0;k_block<uinshape[2];k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uinshape[2]);k++){
#else
                    for (int_t k = 0;k<uinshape[2];k++){
#endif
                            outdev[idx(2*i,2*j,2*k,uouthalo,uoutstride)] = indev[idx(i,j,k,uinhalo,uinstride)];
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(2:2:end,:,:) = (uf(3:2:end,:,:) + uf(1:2:end-2,:,:))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 1;i<uoutshape[0]-1;i+=2){
                for (int_t j = 0;j<uoutshape[1];j+=2){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 0;k_block<uoutshape[2];k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uoutshape[2]);k+=2){
#else
                    for (int_t k = 0;k<uoutshape[2];k+=2){
#endif
                            outdev[idx(i,j,k,uouthalo,uoutstride)] = 0.5*(outdev[idx(i-1,j,k,uouthalo,uoutstride)]+outdev[idx(i+1,j,k,uouthalo,uoutstride)]);
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(:,2:2:end,:) = (uf(:,3:2:end,:) + uf(:,1:2:end-2,:))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uoutshape[0];i++){
                for (int_t j = 1;j<uoutshape[1]-1;j+=2){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 0;k_block<uoutshape[2];k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uoutshape[2]);k+=2){
#else
                    for (int_t k = 0;k<uoutshape[2];k+=2){
#endif
                            outdev[idx(i,j,k,uouthalo,uoutstride)] = 0.5*(outdev[idx(i,j-1,k,uouthalo,uoutstride)]+outdev[idx(i,j+1,k,uouthalo,uoutstride)]);
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
        #pragma omp target device(uout.device) is_device_ptr(outdev)
        {
            //uf(:,:,2:2:end) = (uf(:,:,3:2:end) + uf(:,:,1:2:end-2))/2;
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = 0;i<uoutshape[0];i++){
                for (int_t j = 0;j<uoutshape[1];j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = 1;k_block<uoutshape[2]-1;k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,uoutshape[2]-1);k+=2){
#else
                    for (int_t k = 1;k<uoutshape[2]-1;k+=2){
#endif
                            outdev[idx(i,j,k,uouthalo,uoutstride)] = 0.5*(outdev[idx(i,j,k-1,uouthalo,uoutstride)]+outdev[idx(i,j,k+1,uouthalo,uoutstride)]);
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
