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
        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,omega)
        {
            // Red Points
            #pragma omp teams distribute parallel for collapse(3) schedule(static,8)
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
                    for (int_t k = zmin;k<zmax;k++){
                        if (u.idx(i,j,k) % 2 == 0)
                        udev[u.idx(i,j,k)] = (1.0-omega)*udev[u.idx(i,j,k)] + (omega/6.0)*(udev[u.idx((i-1),j,k)] + udev[u.idx((i+1),j,k)]
                                                +udev[u.idx(i,(j-1),k)] + udev[u.idx(i,(j+1),k)]
                                                +udev[u.idx(i,j,(k-1))] + udev[u.idx(i,j,(k+1))] - hsq*fdev[f.idx(i,j,k)]);
                    }
                }
            }
        }
        #pragma omp target device(u.device) is_device_ptr(udev,fdev) firstprivate(hsq,omega)
        {
            #pragma omp teams distribute parallel for collapse(3) schedule(static,8)
            for (int_t i = xmin;i<xmax;i++){
                for (int_t j = ymin;j<ymax;j++){
                    for (int_t k = zmin;k<zmax;k++){
                        if (u.idx(i,j,k) % 2 == 1)
                        udev[u.idx(i,j,k)] = (1.0-omega)*udev[u.idx(i,j,k)] + (omega/6.0)*(udev[u.idx((i-1),j,k)] + udev[u.idx((i+1),j,k)]
                                                +udev[u.idx(i,(j-1),k)] + udev[u.idx(i,(j+1),k)]
                                                +udev[u.idx(i,j,(k-1))] + udev[u.idx(i,j,(k+1))] - hsq*fdev[f.idx(i,j,k)]);
                    }
                }
            }
        }
    }
}
#endif
