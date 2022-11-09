#import "array.h"
#import "boundary.h"
#import <iostream>
#import <cstdlib>
#import <cassert>

using std::cout;
using std::endl;

int main(void){
    cout << "Devices available: " << omp_get_num_devices() << endl;
    DeviceArray<double_t>A(omp_get_default_device(),{10,2,1});
    for (int i=0;i<A.size;i++){
        A.at[i] = i;
    }
    A.to_device();
    double_t * Ad = A.devptr;
    #pragma omp target teams distribute parallel for schedule(static) is_device_ptr(Ad) device(A.device)
    for (int i=0;i<A.size;i++){
        Ad[i] = 2*Ad[i];
    }
    A.to_host();
    for (int i=0;i<A.size;i++){
        assert(A.at[i]==2*i);
    }
    cout << "DeviceArray works as intended" << endl;

    return EXIT_SUCCESS;
}