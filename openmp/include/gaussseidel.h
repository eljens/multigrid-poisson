#ifndef POISSON_GAUSS_SEIDEL
#define POISSON_GAUSS_SEIDEL

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"

#include <iostream>

namespace Poisson{
    template <typename T>
    class Domain;

    template <class T>
    class GaussSeidel : 
        public Relaxation<T> {
            private:
                bool alternate = true;
            public:
                void relax(Domain<T>& domain,T omega);
                constexpr T default_omega();
                constexpr int smoothing_multiplier() {return 2;};
                constexpr bool requires_duplicate_solution();
                void relaxation_kernel(Domain<T>& domain,T omega,const int_t xmin,const int_t xmax,const int_t ymin,const int_t ymax,const int_t zmin,const int_t zmax);
    };

    template <class T>
    void GaussSeidel<T>::relaxation_kernel(Domain<T>& domain,T omega,const int_t xmin,const int_t xmax,const int_t ymin,const int_t ymax,const int_t zmin,const int_t zmax){
        DeviceArray<T>& u = *domain.u;
        const DeviceArray<T>& f = *domain.f;
        const T one_omega = 1.0-omega;
        const T omega_sixth = (omega/6.0);
        const T hsq = domain.settings.h*domain.settings.h;
        T * udev = u.devptr;
        T * fdev = f.devptr;

        const Halo & uhalo = u.halo;
        const uint_t (&ustride)[3] = u.stride;

        const Halo & fhalo = f.halo;
        const uint_t (&fstride)[3] = f.stride;
        if (alternate) {
        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,one_omega,omega_sixth)
        {
            // Red Points
            #pragma omp teams distribute parallel for collapse(3) SCHEDULE DIST_SCHEDULE
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = zmin;k_block<zmax;k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,zmax);k++){
#else
                    for (int_t k = zmin;k<zmax;k++){
#endif
                            if (idx(i,j,k,uhalo,ustride) % 2 == 0)
                            udev[idx(i,j,k,uhalo,ustride)] = one_omega*udev[idx(i,j,k,uhalo,ustride)] + omega_sixth*(udev[idx((i-1),j,k,uhalo,ustride)] + udev[idx((i+1),j,k,uhalo,ustride)]
                                                    +udev[idx(i,(j-1),k,uhalo,ustride)] + udev[idx(i,(j+1),k,uhalo,ustride)]
                                                    +udev[idx(i,j,(k-1),uhalo,ustride)] + udev[idx(i,j,(k+1),uhalo,ustride)] - hsq*fdev[idx(i,j,k,fhalo,fstride)]);
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
        }
        else {
        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,one_omega,omega_sixth)
        {
            #pragma omp teams distribute parallel for collapse(3) SCHEDULE DIST_SCHEDULE
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
#ifdef BLOCK_SIZE
                    for (int_t k_block = zmin;k_block<zmax;k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,zmax);k++){
#else
                    for (int_t k = zmin;k<zmax;k++){
#endif
                            if (idx(i,j,k,uhalo,ustride) % 2 == 1)
                            udev[idx(i,j,k,uhalo,ustride)] = one_omega*udev[idx(i,j,k,uhalo,ustride)] + omega_sixth*(udev[idx((i-1),j,k,uhalo,ustride)] + udev[idx((i+1),j,k,uhalo,ustride)]
                                                    +udev[idx(i,(j-1),k,uhalo,ustride)] + udev[idx(i,(j+1),k,uhalo,ustride)]
                                                    +udev[idx(i,j,(k-1),uhalo,ustride)] + udev[idx(i,j,(k+1),uhalo,ustride)] - hsq*fdev[idx(i,j,k,fhalo,fstride)]);
#ifdef BLOCK_SIZE
                        }
#endif
                    }
                }
            }
        }
        }
        alternate = !alternate;
    }

    template <class T>
    constexpr T GaussSeidel<T>::default_omega(){
        return 8.0/9.0;
    }

    template <class T>
    constexpr bool GaussSeidel<T>::requires_duplicate_solution(){
        return false;
    }
}
#endif
