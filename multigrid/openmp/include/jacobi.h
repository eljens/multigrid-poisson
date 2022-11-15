#ifndef POISSON_JACOBI
#define POISSON_JACOBI

#define CHUNK_SIZE 8

#include "domain.h"
#include "devicearray.h"

#include <iostream>

template <typename T>
class Domain;

using std::cout;
using std::endl;

template <class T>
void jacobi(Domain<T>& domain,T omega){
    DeviceArray<T>& u = (domain.even) ? *domain.u : *domain.uprev;
    const DeviceArray<T>& v = (domain.even) ? *domain.uprev : *domain.u;
    const DeviceArray<T>& f = *domain.f;

    const T hsq = domain.settings.h*domain.settings.h;

    T * udev = u.devptr;
    T * vdev = v.devptr;
    T * fdev = f.devptr; 
    #pragma omp target map(hsq,omega) device(u.device) is_device_ptr(udev,vdev,fdev)
    {
        #pragma omp teams distribute parallel for collapse(3) schedule(static,1)
        for (int_t i = 1;i<u.shape[0]-1;i++){
            for (int_t j = 1;j<u.shape[1]-1;j++){
                for (int_t k = 1;k<u.shape[2]-1;k++){
                    udev[u.idx(i,j,k)] = (1.0-omega)*vdev[v.idx(i,j,k)];
                    udev[u.idx(i,j,k)] += (omega/6.0)*(vdev[v.idx((i-1),j,k)] + vdev[v.idx((i+1),j,k)]
                                            +vdev[v.idx(i,(j-1),k)] + vdev[v.idx(i,(j+1),k)]
                                            +vdev[v.idx(i,j,(k-1))] + vdev[v.idx(i,j,(k+1))] - hsq*fdev[f.idx(i,j,k)]);
                }
            }
        }
    }
}

#endif