#ifndef POISSON_GAUSS_SEIDEL
#define POISSON_GAUSS_SEIDEL

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"

#include <iostream>

namespace Poisson{
    template <typename T>
    class Domain;

    using std::cout;
    using std::endl;

    template <class T>
    void gaussseidel(Domain<T>& domain,T omega){

        DeviceArray<T>& u = *domain.u;
        const DeviceArray<T>& f = *domain.f;

        const T hsq = domain.settings.h*domain.settings.h;

        // Updating boundaries
        domain.east->update(u,domain.settings);
        domain.west->update(u,domain.settings);
        domain.north->update(u,domain.settings);
        domain.south->update(u,domain.settings);
        domain.top->update(u,domain.settings);
        domain.bottom->update(u,domain.settings);

        // Red Black Gauss Seidel iteration 
        T * udev = u.devptr;
        T * fdev = f.devptr;

        const int xmin = 1-domain.halo.west;
        const int xmax = u.shape[0]-1+domain.halo.east;
        const int ymin = 1-domain.halo.south;
        const int ymax = u.shape[1]-1+domain.halo.north;
        const int zmin = 1-domain.halo.bottom;
        const int zmax = u.shape[2]-1+domain.halo.top;

        const Halo & uhalo = u.halo;
        const uint_t (&ustride)[3] = u.stride;

        const Halo & fhalo = f.halo;
        const uint_t (&fstride)[3] = f.stride;

        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,omega)
        {
            // Red Points
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
                    for (int_t k = zmin;k<zmax;k++){
                        if (idx(i,j,k,uhalo,ustride) % 2 == 0)
                        udev[idx(i,j,k,uhalo,ustride)] = (1.0-omega)*udev[idx(i,j,k,uhalo,ustride)] + (omega/6.0)*(udev[idx((i-1),j,k,uhalo,ustride)] + udev[idx((i+1),j,k,uhalo,ustride)]
                                                +udev[idx(i,(j-1),k,uhalo,ustride)] + udev[idx(i,(j+1),k,uhalo,ustride)]
                                                +udev[idx(i,j,(k-1),uhalo,ustride)] + udev[idx(i,j,(k+1),uhalo,ustride)] - hsq*fdev[idx(i,j,k,fhalo,fstride)]);
                    }
                }
            }
        }
        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,omega)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
                    for (int_t k = zmin;k<zmax;k++){
                        if (idx(i,j,k,uhalo,ustride) % 2 == 1)
                        udev[idx(i,j,k,uhalo,ustride)] = (1.0-omega)*udev[idx(i,j,k,uhalo,ustride)] + (omega/6.0)*(udev[idx((i-1),j,k,uhalo,ustride)] + udev[idx((i+1),j,k,uhalo,ustride)]
                                                +udev[idx(i,(j-1),k,uhalo,ustride)] + udev[idx(i,(j+1),k,uhalo,ustride)]
                                                +udev[idx(i,j,(k-1),uhalo,ustride)] + udev[idx(i,j,(k+1),uhalo,ustride)] - hsq*fdev[idx(i,j,k,fhalo,fstride)]);
                    }
                }
            }
        }
    }
}
#endif