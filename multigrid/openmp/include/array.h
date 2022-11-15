#ifndef HOST_ARRAY
#define HOST_ARRAY

#include "definitions.h"
#include <iostream>
#include "omp.h"
#include "halo.h"
#include "parser.h"

#include <iostream>
#include <stdio.h>
#include <inttypes.h>

using std::initializer_list;
using std::cout;
using std::cerr;
using std::endl;
using std::string;

template <class T>
class Array {
	private:
		void allocator();
	protected:
	public:
		T * at;
		const uint_t ndims = 3;
		const uint_t size;
		const uint_t shape[3];
		const uint_t stride[3];
		const Halo halo;
		Array(uint_t i, uint_t j, uint_t k);

		Array(uint_t i, uint_t j, uint_t k, Halo & _halo);
	
		virtual ~Array();

		inline uint_t idx(uint_t i,uint_t j, uint_t k) const;

		inline uint_t idx_halo(uint_t i,uint_t j, uint_t k) const;

		void print(Settings & settings,const char * str) const;
};

template <class T>
Array<T>::Array(uint_t i, uint_t j, uint_t k) : size(i*j*k),shape{i,j,k}, stride{j*k,k,1}, halo(0,0,0,0,0,0) {
	this->allocator();
}

template <class T>
Array<T>::Array(uint_t i, uint_t j, uint_t k, Halo & _halo) :
	size((i+_halo.east+_halo.west)*(j+_halo.north+_halo.south)*(k+_halo.top+_halo.bottom)),
	shape{i,j,k},
	stride{(j+_halo.north+_halo.south)*(k+_halo.top+_halo.bottom),(k+_halo.top+_halo.bottom),1}, 
	halo(_halo){
	this->allocator();
}

template <class T>
Array<T>::~Array() {
	delete[] this->at;
}

#pragma omp declare target
template <class T>
inline uint_t Array<T>::idx(uint_t i,uint_t j, uint_t k) const{
	uint_t res = (i+this->halo.west)*this->stride[0] + (j+this->halo.south)*this->stride[1] + (k+this->halo.bottom)*this->stride[2];
	return res;
}

template <class T>
inline uint_t Array<T>::idx_halo(uint_t i,uint_t j, uint_t k) const {
	uint_t res = i*this->stride[0] + j*this->stride[1] + k*this->stride[2];
	return res;
}
#pragma omp end declare target

template <class T>
void Array<T>::allocator(){
	try {
		this->at = new T[size];
	}
	catch (std::bad_alloc&) {
		cerr << "Memory allocation failed in Array" << endl;
	}
}

int is_little_endian(void) {
    int num = 1;
    return (*((char *)&num) == 1);
}

template <class T>
void Array<T>::print(Settings & settings, const char *fname) const{

    FILE *f_ptr;
    uint_t written = 0;
    uint_t items = this->size;
    uint_t i,j,k;

    if ( (f_ptr = fopen(fname, "w")) == NULL ) {
       perror("No output! fopen()");
       return;
    }

    // Write VTK file header
    fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
    fprintf(f_ptr, "saved from function print_vtk.\n");
    fprintf(f_ptr, "BINARY\n");
    fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(f_ptr, "DIMENSIONS %d %d %d\n",(int) this->shape[0],(int) this->shape[1],(int) this->shape[2]);
    fprintf(f_ptr, "ORIGIN %f %f %f\n",(double) settings.origin[0],(double) settings.origin[1],(double) settings.origin[2]);
    fprintf(f_ptr, "SPACING %f %f %f\n",(double) settings.h,(double) settings.h,(double) settings.h);
    fprintf(f_ptr, "POINT_DATA %lu\n", items);
    fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
    fprintf(f_ptr, "LOOKUP_TABLE default\n");

    if ( is_little_endian() ) {
        // System is little endian, so we need to reverse the byte order.
        for (k = 0; k < this->shape[2]; ++k) {
            for (j = 0; j < this->shape[1]; ++j) {
                for (i = 0; i < this->shape[0]; ++i) {
            uint64_t crnt = *(uint64_t *)(&(this->at[this->idx(i,j,k)])); // Get double as int

            // Reverse byte order and write to file
            crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
            crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
            crnt = (crnt & 0x00FF00FF00FF00FF) << 8  | (crnt & 0xFF00FF00FF00FF00) >> 8;
            written += fwrite(&crnt, sizeof(uint64_t), 1, f_ptr);
                }
            }
        }
    } else {
        for (k = 0; k < this->shape[2]; ++k) {
            for (j = 0; j < this->shape[1]; ++j) {
                for (i = 0; i < this->shape[0]; ++i) {
                    written += fwrite(&(this->at[this->idx(i,j,k)]), sizeof(T), 1, f_ptr);
                }
            }
        }
    }

    if ( written != items ) {
	    cerr << "Failed to print file " << fname << endl;
    }

    fclose(f_ptr);
}

#endif
