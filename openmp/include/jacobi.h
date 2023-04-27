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
                void relax(Domain<T>& domain,T omega);
                constexpr T default_omega();
                constexpr bool requires_duplicate_solution();
    };

    template <class T>
    void Jacobi<T>::relax(Domain<T>& domain,T omega){
        domain.swap_u();

        DeviceArray<T>& u = *domain.u;
        DeviceArray<T>& v = *domain.uprev;
        const DeviceArray<T>& f = *domain.f;

        const T hsq = domain.settings.h*domain.settings.h;

        // Updating boundaries
        domain.east->update(v,domain.settings);
        domain.west->update(v,domain.settings);
        domain.north->update(v,domain.settings);
        domain.south->update(v,domain.settings);
        domain.top->update(v,domain.settings);
        domain.bottom->update(v,domain.settings);

        // Jacobi iteration 
        T * udev = u.devptr;
        T * vdev = v.devptr;
        T * fdev = f.devptr;

        const int xmin = 1-domain.halo.west;
        const int xmax = u.shape[0]-1+domain.halo.east;
        const int ymin = 1-domain.halo.south;
        const int ymax = u.shape[1]-1+domain.halo.north;
        const int zmin = 1-domain.halo.bottom;
        const int zmax = u.shape[2]-1+domain.halo.top;

        const Halo & uhalo = u.halo;
        const uint_t (&ustride)[3] = u.stride;

        const Halo & vhalo = v.halo;
        const uint_t (&vstride)[3] = v.stride;

        const Halo & fhalo = f.halo;
        const uint_t (&fstride)[3] = f.stride;

        const T one_omega = 1.0-omega;
        const T omega_sixth = (omega/6.0);

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

        if (domain.east->is_internal_boundary()){
            OMPBoundary<T> * omp_east = (OMPBoundary<T> *) domain.east;
            omp_east->fill_send_buffer(u,domain.settings);
        }

        if (domain.west->is_internal_boundary()){
            OMPBoundary<T> * omp_west = (OMPBoundary<T> *) domain.west;
            omp_west->fill_send_buffer(u,domain.settings);
        }

        //domain.swap_u();
    }

    template <class T>
    constexpr T Jacobi<T>::default_omega(){
        return 6.0/7.0;
    }

    template <class T>
    constexpr bool Jacobi<T>::requires_duplicate_solution(){
        return true;
    }
}
#endif
