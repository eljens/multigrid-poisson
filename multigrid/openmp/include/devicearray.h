#ifndef DEVICE_ARRAY
#define DEVICE_ARRAY

#include "array.h"
#include "parser.h"
#include "halo.h"

template <class T>
class DeviceArray :
	public Array<T> {
		private:
			void device_allocator();
		protected:
			uint host = omp_get_initial_device();
		public:
			uint_t device;
			T * devptr;

			DeviceArray(uint_t _device,uint_t i, uint_t j, uint_t k);

			DeviceArray(Settings & settings, Halo & _halo);

			virtual ~DeviceArray();

			void to_device();

			void to_host();
};

template <class T>
DeviceArray<T>::DeviceArray(uint_t _device,uint_t i, uint_t j, uint_t k) : 
	Array<T>(i,j,k), device(_device) {
	this->device_allocator();
}

template <class T>
DeviceArray<T>::DeviceArray(Settings & settings, Halo & _halo) :
	Array<T>(
		settings.dims[0],
		settings.dims[1],
		settings.dims[2],
		_halo), device(settings.dev) {
	this->device_allocator();
}

template <class T>
DeviceArray<T>::~DeviceArray(){
	omp_target_free(this->devptr,this->device);
	#pragma omp target exit data map(delete:this->stride[:this->ndims],this->shape[:this->ndims],this->halo) device(this->device)
}

template <class T>
void DeviceArray<T>::to_device() {
	int_t res = omp_target_memcpy(this->devptr,this->at,(int) this->size*sizeof(T),0,0,(int) this->device,(int) this->host);
	if (res != 0){
		cerr << "Error on device " << this->device << ": omp_target_memcpy returned " << res << endl;
	}
}

template <class T>
void DeviceArray<T>::to_host() {
	int_t res = omp_target_memcpy(this->at,this->devptr,(int) this->size*sizeof(T),0,0,(int) this->host,(int) this->device);
	if (res != 0){
		cerr <<  "Error on device " << this->device << ":omp_target_memcpy returned " << res << endl;
	}
}

template <class T>
void DeviceArray<T>::device_allocator() {
	this->devptr = (T*) omp_target_alloc((int)this->size*sizeof(T),(int)this->device);
	if (this->devptr == NULL){
		cerr << "Error allocating DeviceArray on device " << this->device << endl;
	}
	#pragma omp target enter data map(to:this->stride[:this->ndims],this->shape[:this->ndims],this->halo) device(this->device)		
}

#endif
