#ifndef POISSON_JACOBI
#define POISSON_JACOBI

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"
#include "relaxation.h"

#include <iostream>

namespace Poisson{
    template <typename T>
    class Domain;

    using std::cout;
    using std::endl;

    template <class T>
    class Jacobi : 
        public Relaxation<T> {
            public:
                constexpr T default_omega();
                constexpr bool requires_duplicate_solution();
                constexpr int smoothing_multiplier() {return 1;};
                void relaxation_kernel(Domain<T>& domain,T omega,const int_t xmin,const int_t xmax,const int_t ymin,const int_t ymax,const int_t zmin,const int_t zmax);
    };

    template <class T>
    constexpr T Jacobi<T>::default_omega(){
        return 6.0/7.0;
    }

    template <class T>
    constexpr bool Jacobi<T>::requires_duplicate_solution(){
        return true;
    }

    template <class T>
    void Jacobi<T>::relaxation_kernel(Domain<T>& domain,T omega,const int_t xmin,const int_t xmax,const int_t ymin,const int_t ymax,const int_t zmin,const int_t zmax){
        DeviceArray<T>& u = *domain.u;
        DeviceArray<T>& v = *domain.uprev;
        const DeviceArray<T>& f = *domain.f;
        const T one_omega = 1.0-omega;
        const T omega_sixth = (omega/6.0);
        const T hsq = domain.settings.h*domain.settings.h;
        T * udev = u.devptr;
        T * vdev = v.devptr;
        T * fdev = f.devptr;

        const Halo & uhalo = u.halo;
        const uint_t (&ustride)[3] = u.stride;

        const Halo & vhalo = v.halo;
        const uint_t (&vstride)[3] = v.stride;

        const Halo & fhalo = f.halo;
        const uint_t (&fstride)[3] = f.stride;

        #pragma omp target device(u.device) is_device_ptr(udev,vdev,fdev)
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
                            udev[idx(i,j,k,uhalo,ustride)] = one_omega*vdev[idx(i,j,k,vhalo,vstride)];
                            udev[idx(i,j,k,uhalo,ustride)] += omega_sixth*(vdev[idx((i-1),j,k,vhalo,vstride)] + vdev[idx((i+1),j,k,vhalo,vstride)]
                                                    +vdev[idx(i,(j-1),k,vhalo,vstride)] + vdev[idx(i,(j+1),k,vhalo,vstride)]
                                                    +vdev[idx(i,j,(k-1),vhalo,vstride)] + vdev[idx(i,j,(k+1),vhalo,vstride)] - hsq*fdev[idx(i,j,k,fhalo,fstride)]);
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
