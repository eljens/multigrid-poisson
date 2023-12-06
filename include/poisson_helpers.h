#ifndef MG_POISSON_HELPERS
#define MG_POISSON_HELPERS
#include "libpoisson.h"

using Poisson::DeviceArray;
using Poisson::Boundary;

template<class T>
void update_rhs(DeviceArray<T> & f,const int sizex, const int sizey, const int sizez, const double dx, double* dens_i, double* dens_e, double eps){
	double dxdx = dx * dx;
    double e_0 = 1.602176487E-19;
    #pragma omp parallel for collapse(3) shared(f)
	for(int i = 0; i < sizex; i++){
        for(int j = 0; j < sizey; j++){
            for(int k = 0; k < sizez; k++){
                const int idx = i * (sizey * sizez) + j * sizez + k;
                //note that the Poisson equation to be solved is \Delta \Phi = -f, hence f needs to be rho/eps_0
                //dens[.] is the number density, both with a positive sign
                f.at[f.idx(i,j,k)] = - e_0 * ((dens_i[idx]-dens_e[idx]) / eps);//*dxdx;
            }
        }
    }
}

template<class T>
void write_dirichlet_conditions(Boundary<T> & top,Boundary<T> & bottom, const int sizex, const int sizey, const int sizez, const double dx){
    #pragma omp parallel for collapse(2) shared(background_field,top,bottom)
	for(int i = 0; i < sizex; i++){
        for(int j = 0; j < sizey; j++){
            bottom.arr.at[bottom.arr.idx(i,j,0)] = 0.0;
            top.arr.at[top.arr.idx(i,j,0)] = background_field[2] * dx * sizez;
        }
    }
}

#endif
