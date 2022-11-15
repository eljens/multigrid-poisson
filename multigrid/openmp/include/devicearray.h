#ifndef DEVICE_ARRAY
#define DEVICE_ARRAY

#include "array.h"

template <class T>
class DeviceArray :
	public Array<T> {
		protected:
			uint host = omp_get_initial_device();
		public:
			uint_t device;
			T * devptr;

			DeviceArray(uint_t device,uint_t i, uint_t j, uint_t k);

			virtual ~DeviceArray();

			void to_device();

			void to_host();
};

template <class T>
DeviceArray<T>::DeviceArray(uint_t device,uint_t i, uint_t j, uint_t k) : Array<T>(i,j,k){
	this->device = device;
	this->devptr = (T*) omp_target_alloc((int)this->size*sizeof(T),(int)this->device);
	if (this->devptr == NULL){
		cerr << "Error allocating DeviceArray on device " << this->device << endl;
	}
	#pragma omp target enter data map(to:this->stride[:this->ndims],this->shape[:this->ndims]) device(this->device)
}

template <class T>
DeviceArray<T>::~DeviceArray(){
	omp_target_free(this->devptr,this->device);
	#pragma omp target exit data map(delete:this->stride[:this->ndims],this->shape[:this->ndims]) device(this->device)
}

template <class T>
void DeviceArray<T>::to_device() {
	int_t res = omp_target_memcpy(this->devptr,this->at,(int) this->size*sizeof(T),0,0,(int) this->device,(int) this->host);
	if (res != 0){
		cerr << "omp_target_memcpy returned " << res << endl;
	}
}

template <class T>
void DeviceArray<T>::to_host() {
	int_t res = omp_target_memcpy(this->at,this->devptr,(int) this->size*sizeof(T),0,0,(int) this->host,(int) this->device);
	if (res != 0){
		cerr << "omp_target_memcpy returned " << res << endl;
	}
}

#endif
