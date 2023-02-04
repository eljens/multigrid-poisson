#ifndef POISSON_RESIDUAL
#define POISSON_RESIDUAL

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"

#include <iostream>

using std::cout;
using std::endl;

namespace Poisson{
    template <typename T>
    class Domain;

    template <class T>
    void residual(Domain<T>& domain){
        DeviceArray<T>& u = *domain.u;
        DeviceArray<T>& r = *domain.r;
        const DeviceArray<T>& f = *domain.f;

        const T hsq = domain.settings.h*domain.settings.h;

        // Updating boundaries
        domain.east->update(u,domain.settings);
        domain.west->update(u,domain.settings);
        domain.north->update(u,domain.settings);
        domain.south->update(u,domain.settings);
        domain.top->update(u,domain.settings);
        domain.bottom->update(u,domain.settings);

        // Extracting device pointers
        T * udev = u.devptr;
        T * rdev = r.devptr;
        T * fdev = f.devptr;

        const int xmin = 1-domain.halo.west;
        const int xmax = r.shape[0]-1+domain.halo.east;
        const int ymin = 1-domain.halo.south;
        const int ymax = r.shape[1]-1+domain.halo.north;
        const int zmin = 1-domain.halo.bottom;
        const int zmax = r.shape[2]-1+domain.halo.top;

        const Halo & rhalo = r.halo;
        const uint_t (&rstride)[3] = r.stride;

        const Halo & fhalo = f.halo;
        const uint_t (&fstride)[3] = f.stride;

        const Halo & uhalo = u.halo;
        const uint_t (&ustride)[3] = u.stride;

        #pragma omp target device(u.device) is_device_ptr(udev,rdev,fdev) firstprivate(hsq)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
                    for (int_t k_block = zmin;k_block<zmax;k_block+=BLOCK_SIZE){
                        #pragma omp simd
                        for (int_t k = k_block;k<MIN(k_block+BLOCK_SIZE,zmax);k++){
                            rdev[idx(i,j,k,rhalo,rstride)] = fdev[idx(i,j,k,fhalo,fstride)] - (udev[idx((i-1),j,k,uhalo,ustride)] + udev[idx((i+1),j,k,uhalo,ustride)]
                                                +udev[idx(i,(j-1),k,uhalo,ustride)] + udev[idx(i,(j+1),k,uhalo,ustride)]
                                                +udev[idx(i,j,(k-1),uhalo,ustride)] + udev[idx(i,j,(k+1),uhalo,ustride)] - 6.0*udev[idx(i,j,k,uhalo,ustride)])/hsq;
                        }
                    }
                }
            }
        }
    }
}
#endif
