#include "../include/devicearray.h"
#include <iostream>
#include <cstdlib>
#include <cassert>

using std::cout;
using std::endl;

int main(void){
    cout << "Devices available: " << omp_get_num_devices() << endl;
    DeviceArray<double_t>A(omp_get_default_device(),10,3,2);
    for (int i=0;i<10;i++){
        for (int j=0;j<3;j++){
            for (int k=0;k<2;k++){
                A.at[A.idx(i,j,k)] = i*3*2 + j*2 + k;
            }
        }
    }
    A.to_device();
    double_t * Ad = A.devptr;
    #pragma omp target teams distribute parallel for schedule(static) is_device_ptr(Ad) device(A.device)
    for (int i=0;i<10;i++){
        for (int j=0;j<3;j++){
            for (int k=0;k<2;k++){
                Ad[A.idx(i,j,k)] = 2*Ad[A.idx(i,j,k)];
            }
        }
    }
    A.to_host();
    for (int i=0;i<A.size;i++){
        assert(A.at[i]==2*i);
    }
    cout << "DeviceArray works as intended" << endl;

    return EXIT_SUCCESS;
}