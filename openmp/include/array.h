#ifndef HOST_ARRAY
#define HOST_ARRAY

#include "definitions.h"
#include <iostream>
#include "omp.h"
#include "halo.h"
#include "settings.h"

#include <iostream>
#include <stdio.h>
#include <inttypes.h>

using std::initializer_list;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
namespace Poisson{
    constexpr uint_t idx(const uint_t i,const uint_t j,const uint_t k,const Halo & _halo,const uint_t (&_stride)[3]) 
    {
        uint_t res = (i+_halo.west)*_stride[0] + (j+_halo.south)*_stride[1] + (k+_halo.bottom)*_stride[2];
        return res;
    }

    constexpr uint_t idx_halo(const uint_t i,const uint_t j, const uint_t k, const uint_t (&_stride)[3]) 
    {
        uint_t res = i*_stride[0] + j*_stride[1] + k*_stride[2];
        return res;
    }

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

            constexpr uint_t idx(const uint_t i,const uint_t j,const uint_t k) const;

            constexpr uint_t idx_halo(const uint_t i,const uint_t j, const uint_t k) const;

            void init_zero();

            void add(const Array<T> & arr);

            T infinity_norm() const;

            void print(Settings & settings,const char * str) const;

            void print_halo(Settings & settings,const char * str) const;
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

    template <class T>
    constexpr uint_t Array<T>::idx(const uint_t i,const uint_t j,const uint_t k) const{
        return Poisson::idx(i,j,k,this->halo,this->stride);
    }


    template <class T>
    constexpr uint_t Array<T>::idx_halo(const uint_t i, const uint_t j, const uint_t k) const {
        return Poisson::idx_halo(i,j,k,this->stride);
    }

    template <class T>
    void Array<T>::allocator(){
        try {
            this->at = new T[size];
        }
        catch (std::bad_alloc&) {
            cerr << "Memory allocation failed in Array" << endl;
        }
    }

    template<class T>
    void Array<T>::init_zero(){
        #pragma omp parallel for schedule(static,CHUNK_SIZE)
        for(uint_t i = 0;i<this->size;i++){
            this->at[i] = 0.0;
        }
    }

    template<class T>
    void Array<T>::add(const Array<T> & arr){
        if (!(arr.size == this->size)){
            cerr << "Array size does not match in Array<T>.add()" << endl;
        }
        #pragma omp parallel for collapse(3) schedule(static,CHUNK_SIZE)
        for(uint_t i = 0;i<this->shape[0];i++){
            for(uint_t j = 0;j<this->shape[1];j++){
                for(uint_t k = 0;k<this->shape[2];k++){
                    this->at[this->idx(i,j,k)] += arr.at[arr.idx(i,j,k)];
                }
            }
        }
    }

    template<class T>
    T Array<T>::infinity_norm() const{
        T res = 0.0;
        #pragma omp parallel for collapse(3) schedule(static,CHUNK_SIZE) reduction(max:res)
        for(uint_t i = 0;i<this->shape[0];i++){
            for(uint_t j = 0;j<this->shape[1];j++){
                for(uint_t k = 0;k<this->shape[2];k++){
                    T abselem = std::abs(this->at[this->idx(i,j,k)]);
                    res = std::max(res,abselem);
                }
            }
        }
        return res;
    }

    int is_little_endian(void) {
        int num = 1;
        return (*((char *)&num) == 1);
    }

    template <class T>
    void Array<T>::print(Settings & settings, const char *fname) const{

        FILE *f_ptr;
        uint_t written = 0;
        uint_t items = this->shape[0]*this->shape[1]*this->shape[2];
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

    template <class T>
    void Array<T>::print_halo(Settings & settings, const char *fname) const{

        FILE *f_ptr;
        uint_t written = 0;
        uint_t items = this->size;
        uint_t i,j,k;

        if ( (f_ptr = fopen(fname, "w")) == NULL ) {
        perror("No output! fopen()");
        return;
        }

        uint_t nx = this->shape[0] + this->halo.east + this->halo.west;
        uint_t ny = this->shape[1] + this->halo.north + this->halo.south;
        uint_t nz = this->shape[2] + this->halo.top + this->halo.bottom;

        cout << "Printing file " << fname << " with shape (" << nx << "," << ny << "," << nz << ")" << endl;

        // Write VTK file header
        fprintf(f_ptr, "# vtk DataFile Version 3.0\n");
        fprintf(f_ptr, "saved from function print_vtk.\n");
        fprintf(f_ptr, "BINARY\n");
        fprintf(f_ptr, "DATASET STRUCTURED_POINTS\n");
        fprintf(f_ptr, "DIMENSIONS %d %d %d\n",(int) nx,(int) ny,(int) nz);
        fprintf(f_ptr, "ORIGIN %f %f %f\n",
            (double) settings.origin[0]-((double)this->halo.west)*settings.h,
            (double) settings.origin[1]-((double)this->halo.south)*settings.h,
            (double) settings.origin[2]-((double)this->halo.bottom)*settings.h);
        fprintf(f_ptr, "SPACING %f %f %f\n",(double) settings.h,(double) settings.h,(double) settings.h);
        fprintf(f_ptr, "POINT_DATA %lu\n", items);
        fprintf(f_ptr, "SCALARS %s %s 1\n", "gray", "double");
        fprintf(f_ptr, "LOOKUP_TABLE default\n");

        if ( is_little_endian() ) {
            // System is little endian, so we need to reverse the byte order.
            for (k = 0; k < nz; ++k) {
                for (j = 0; j < ny; ++j) {
                    for (i = 0; i < nx; ++i) {
                uint64_t crnt = *(uint64_t *)(&(this->at[this->idx_halo(i,j,k)])); // Get double as int

                // Reverse byte order and write to file
                crnt = (crnt & 0x00000000FFFFFFFF) << 32 | (crnt & 0xFFFFFFFF00000000) >> 32;
                crnt = (crnt & 0x0000FFFF0000FFFF) << 16 | (crnt & 0xFFFF0000FFFF0000) >> 16;
                crnt = (crnt & 0x00FF00FF00FF00FF) << 8  | (crnt & 0xFF00FF00FF00FF00) >> 8;
                written += fwrite(&crnt, sizeof(uint64_t), 1, f_ptr);
                    }
                }
            }
        } else {
            for (k = 0; k < nz; ++k) {
                for (j = 0; j < ny; ++j) {
                    for (i = 0; i < nx; ++i) {
                        written += fwrite(&(this->at[this->idx_halo(i,j,k)]), sizeof(T), 1, f_ptr);
                    }
                }
            }
        }

        if ( written != items ) {
            cerr << "Failed to print file " << fname << endl;
        }

        fclose(f_ptr);
    }
}

#endif
