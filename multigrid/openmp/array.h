#ifndef POISSON_ARRAY
#define POISSON_ARRAY

#include "definitions.h"
#import <iostream>
#import "omp.h"

using std::initializer_list;
using std::cout;
using std::cerr;
using std::endl;

template <class T>
class Array {
	protected:
	public:
		T * at;
		uint_t * shape;
		uint_t ndims = 0;
		uint_t size;
		Array(initializer_list<uint_t> args) {
			// Counting the number of input arguments
			for (auto arg : args){
				ndims++;
			}

			// Extracting the shape of the array
			this->shape = new uint_t[ndims];
			uint_t i = 0;
			uint_t prod = 1;
			for (auto arg : args){
				shape[i] = arg;
				prod *=arg;
				i++;
			}
			this->size = prod;

			// Allocating the space
			this->at = new T[size];
		}
	
		~Array() {
			delete[] this->shape;
			delete[] this->at;
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
			DeviceArray(uint_t device,initializer_list<uint_t> args) : Array<T>(args){
				this->device = device;
				this->devptr = (T*) omp_target_alloc((int)this->size*sizeof(T),(int)this->device);
				if (this->devptr == NULL){
					cerr << "Error allocating DeviceArray on device " << this->device << endl;
				}
			}

			~DeviceArray(){
				omp_target_free(this->devptr,this->device);
			}

			void to_device(){
				int_t res = omp_target_memcpy(this->devptr,this->at,(int) this->size*sizeof(T),0,0,(int) this->device,(int) this->host);
				if (res != 0){
					cerr << "omp_target_memcpy returned " << res << endl;
				}
			}

			void to_host(){
				int_t res = omp_target_memcpy(this->at,this->devptr,(int) this->size*sizeof(T),0,0,(int) this->host,(int) this->device);
				if (res != 0){
					cerr << "omp_target_memcpy returned " << res << endl;
				}
			}
};

#endif

