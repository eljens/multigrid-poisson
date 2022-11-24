#ifndef POISSON_JACOBI
#define POISSON_JACOBI

#include "definitions.h"
#include "domain.h"
#include "devicearray.h"

#include <iostream>

template <typename T>
class Domain;

using std::cout;
using std::endl;

template <class T>
void jacobi(Domain<T>& domain,T omega){
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

    #pragma omp target device(u.device) is_device_ptr(udev,vdev,fdev) firstprivate(hsq,omega)
    {
        #pragma omp teams distribute parallel for collapse(3) schedule(static,8)
        for (int_t i = xmin;i<xmax;i++){
            for (int_t j = ymin;j<ymax;j++){
                for (int_t k = zmin;k<zmax;k++){
                    udev[u.idx(i,j,k)] = (1.0-omega)*vdev[v.idx(i,j,k)];
                    udev[u.idx(i,j,k)] += (omega/6.0)*(vdev[v.idx((i-1),j,k)] + vdev[v.idx((i+1),j,k)]
                                            +vdev[v.idx(i,(j-1),k)] + vdev[v.idx(i,(j+1),k)]
                                            +vdev[v.idx(i,j,(k-1))] + vdev[v.idx(i,j,(k+1))] - hsq*fdev[f.idx(i,j,k)]);
                }
            }
        }
    }
    //domain.swap_u();
}

#endif