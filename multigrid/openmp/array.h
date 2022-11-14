#ifndef POISSON_ARRAY
#define POISSON_ARRAY

#include "definitions.h"
#include <iostream>
#include "omp.h"
extern "C" {
	#include "print.h"
}

#include <iostream>

using std::initializer_list;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

template <class T>
class Array {
	protected:
	public:
		T * at;
		const uint_t shape[3];
		uint_t ndims = 3;
		uint_t size;
		const uint_t stride[3];
		Array(uint_t i, uint_t j, uint_t k) : stride{j*k,k,1},shape{i,j,k} {
			this->size = i*j*k;

			// Allocating the space
			try {
				this->at = new T[size];
			}
			catch (std::bad_alloc&) {
				cerr << "Memory allocation failed in Array" << endl;
			}
		}
	
		virtual ~Array() {
			delete[] this->at;
		}
		#pragma omp declare target
		inline uint_t idx(uint_t i,uint_t j, uint_t k) const{
			uint_t res = i*this->stride[0] + j*this->stride[1] + k*this->stride[2];
			/*if (res >= this->size){
				cerr << "Indexing out of bounds: " << res << " > " << this->size << endl;
			}*/
			return res;
		}
		#pragma omp end declare target

		void to_vtk_file(const char * str) const {
			print_vtk(str,this->shape[0],this->shape[1],this->shape[2],this->at);
			cout << "Printed file " << str << endl;
		}
};

template <class T>
class DeviceArray :
	public Array<T> {
		protected:
			uint host = omp_get_initial_device();
		public:
			uint_t device;
			T * devptr;
			DeviceArray(uint_t device,uint_t i, uint_t j, uint_t k) : Array<T>(i,j,k){
				this->device = device;
				this->devptr = (T*) omp_target_alloc((int)this->size*sizeof(T),(int)this->device);
				if (this->devptr == NULL){
					cerr << "Error allocating DeviceArray on device " << this->device << endl;
				}
				#pragma omp target enter data map(to:this->stride[:this->ndims],this->shape[:this->ndims]) device(this->device)
			}

			virtual ~DeviceArray(){
				omp_target_free(this->devptr,this->device);
			}

			void to_device() {
				int_t res = omp_target_memcpy(this->devptr,this->at,(int) this->size*sizeof(T),0,0,(int) this->device,(int) this->host);
				if (res != 0){
					cerr << "omp_target_memcpy returned " << res << endl;
				}
			}

			void to_host() {
				int_t res = omp_target_memcpy(this->at,this->devptr,(int) this->size*sizeof(T),0,0,(int) this->host,(int) this->device);
				if (res != 0){
					cerr << "omp_target_memcpy returned " << res << endl;
				}
			}
};

#endif

