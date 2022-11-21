#ifndef POISSON_INJECTION
#define POISSON_INJECTION

#include "restriction.h"

template<class T>
class Injection :
    public Restriction<T> {
        public:
        Injection();

        ~Injection();

        void restrict(DeviceArray<T> & uin,DeviceArray<T> & uout);
};

template<class T>
Injection<T>::Injection() : Restriction<T>() {
    // Does nothing
}

template<class T>
Injection<T>::~Injection(){
    // Does nothing
}

template<class T>
void Injection<T>::restrict(DeviceArray<T> & uin,DeviceArray<T> & uout){
    bool cond0 = ((int_t) (uin.shape[0]-1)/2)+1 == uout.shape[0];
    bool cond1 = ((int_t) (uin.shape[1]-1)/2)+1 == uout.shape[1];
    bool cond2 = ((int_t) (uin.shape[2]-1)/2)+1 == uout.shape[2];
    if ((!cond0) || (!cond1) || (!cond2)){
        cerr << "Restriction from domain of size (" << uin.shape[0] << "," << uin.shape[1] << "," << uin.shape[2] << ")";
        cerr << " to (" << uout.shape[0] << "," << uout.shape[1] << "," << uout.shape[2] << ") is incompatible with injection" << endl;
    }

    T * indev = uin.devptr;
    T * outdev = uout.devptr;

    #pragma omp target device(uout.device) is_device_ptr(indev,outdev)
    {
        #pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
        for (int_t i = 0;i<uout.shape[0];i++){
            for (int_t j = 0;j<uout.shape[1];j++){
                for (int_t k = 0;k<uout.shape[2];k++){
                    outdev[uout.idx(i,j,k)] = indev[uin.idx(2*i,2*j,2*k)];
                }
            }
        }
    }
}

#endif