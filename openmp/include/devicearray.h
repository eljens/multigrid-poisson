#ifndef DEVICE_ARRAY
#define DEVICE_ARRAY

#include <exception>

#include "array.h"
#include "settings.h"
#include "halo.h"

namespace Poisson{
	template <class T>
	class DeviceArray :
		public Array<T> {
			private:
				void device_allocator();
				bool on_device = false;
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

				void init_zero();

				void add(const DeviceArray<T> & arr);

				T infinity_norm() const;
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
		this->on_device = true;
	}

	template <class T>
	void DeviceArray<T>::to_host() {
		int_t res = omp_target_memcpy(this->at,this->devptr,(int) this->size*sizeof(T),0,0,(int) this->host,(int) this->device);
		if (res != 0){
			cerr <<  "Error on device " << this->device << ":omp_target_memcpy returned " << res << endl;
		}
		//this->on_device = false;
	}

	template <class T>
	void DeviceArray<T>::device_allocator() {
		if (omp_get_num_devices() < 1){
			throw std::runtime_error("Device array requested but no devices were detected by the OpenMP runtime\n");
		}
		this->devptr = (T*) omp_target_alloc((int)this->size*sizeof(T),(int)this->device);
		if (this->devptr == NULL){
			cerr << "Error allocating DeviceArray on device " << this->device << endl;
		}
		#pragma omp target enter data map(to:this->stride[:this->ndims],this->shape[:this->ndims],this->halo) device(this->device)		
	}

	template<class T>
	void DeviceArray<T>::init_zero(){
		if (!(this->on_device)){
			Array<T>::init_zero();
		}
		else {
			const uint_t _dev = this->device;
			const uint_t _size = this->size;
			T * _devptr = this->devptr;
			#pragma omp target device(_dev) is_device_ptr(_devptr) firstprivate(_size)
			{
				#pragma omp teams distribute parallel for schedule(static,CHUNK_SIZE)
				for(uint_t i = 0;i<_size;i++){
					_devptr[i] = 0.0;
				}
			}
		}
	}

	template<class T>
	void DeviceArray<T>::add(const DeviceArray<T> & arr){
		if (!(this->on_device)){
			Array<T>::add(arr);
		}
		else {
			const uint_t _dev = this->device;
			T * arrdev = arr.devptr;
			T * _devptr = this->devptr;
			const uint_t * _shape = this->shape;
			#pragma omp target device(_dev) is_device_ptr(_devptr,arrdev)
			{
				#pragma omp teams distribute parallel for collapse(3) schedule(static,CHUNK_SIZE)
				for(uint_t i = 0;i<_shape[0];i++){
					for(uint_t j = 0;j<_shape[1];j++){
						for(uint_t k = 0;k<_shape[2];k++){
							_devptr[this->idx(i,j,k)] += arrdev[arr.idx(i,j,k)];
						}
					}
				}
			}
		}
	}

	template<class T>
	T DeviceArray<T>::infinity_norm() const{
		T res = 0.0;
		if (!(this->on_device)){
			res =  Array<T>::infinity_norm();
			cout << "called Array<T>::infinity_norm()" << endl;
		}
		else {
			const uint_t _dev = this->device;
			const uint_t * _shape = this->shape;
			T * _devptr = this->devptr;
			#pragma omp target device(_dev) is_device_ptr(_devptr) map(always,tofrom:res)
			{
				#pragma omp teams distribute parallel for reduction(max:res) collapse(3) schedule(static,CHUNK_SIZE)
				for(uint_t i = 0;i<_shape[0];i++){
					for(uint_t j = 0;j<_shape[1];j++){
						for(uint_t k = 0;k<_shape[2];k++){
							T abselem = std::abs(_devptr[this->idx(i,j,k)]);
							res = std::max(abselem,res);
						}
					}
				}
			}
		}
		return res;
	}
}
#endif
