#ifndef POISSON_RESIDUAL
#define POISSON_RESIDUAL

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"

#include <iostream>

template <typename T>
class Domain;

using std::cout;
using std::endl;

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

    #pragma omp target device(u.device) is_device_ptr(udev,rdev,fdev) firstprivate(hsq)
    {
        #pragma omp teams distribute parallel for collapse(3) schedule(static,8)
        for (int_t i = xmin;i<xmax;i++){
            for (int_t j = ymin;j<ymax;j++){
                for (int_t k = zmin;k<zmax;k++){
                    rdev[r.idx(i,j,k)] = fdev[f.idx(i,j,k)] - (udev[u.idx((i-1),j,k)] + udev[u.idx((i+1),j,k)]
                                            +udev[u.idx(i,(j-1),k)] + udev[u.idx(i,(j+1),k)]
                                            +udev[u.idx(i,j,(k-1))] + udev[u.idx(i,j,(k+1))] - 6.0*udev[u.idx(i,j,k)])/hsq;
                }
            }
        }
    }
}

#endif