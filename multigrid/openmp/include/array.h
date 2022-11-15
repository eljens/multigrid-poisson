#ifndef HOST_ARRAY
#define HOST_ARRAY

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
		const uint_t ndims = 3;
		const uint_t size;
		const uint_t shape[3];
		const uint_t stride[3];
		Array(uint_t i, uint_t j, uint_t k);
	
		virtual ~Array() {
			delete[] this->at;
		}

		inline uint_t idx(uint_t i,uint_t j, uint_t k) const;

		void to_vtk_file(const char * str) const;
};

template <class T>
Array<T>::Array(uint_t i, uint_t j, uint_t k) : size(i*j*k),shape{i,j,k}, stride{j*k,k,1} {
	// Allocating the space
	try {
		this->at = new T[size];
	}
	catch (std::bad_alloc&) {
		cerr << "Memory allocation failed in Array" << endl;
	}
}

#pragma omp declare target
template <class T>
inline uint_t Array<T>::idx(uint_t i,uint_t j, uint_t k) const{
	uint_t res = i*this->stride[0] + j*this->stride[1] + k*this->stride[2];
	return res;
	}
#pragma omp end declare target

template <class T>
void Array<T>::to_vtk_file(const char * str) const {
	print_vtk(str,this->shape[0],this->shape[1],this->shape[2],this->at);
	cout << "Printed file " << str << endl;
}

#endif
