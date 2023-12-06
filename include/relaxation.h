#ifndef POISSON_RELAXATION
#define POISSON_RELAXATION

#include "domain.h"

namespace Poisson{
    template <class T>
    class Relaxation {
            public:
            void relax(Domain<T> *** domain,T omega,int_t nsmooth,int_t level,int_t num_devices);
            virtual T default_omega()=0;
            virtual bool requires_duplicate_solution()=0;
            virtual int smoothing_multiplier()=0;
            void fill_send_buffer(Domain<T>& domain, const T omega, Location_t loc);
            virtual void relaxation_kernel(Domain<T>& domain,T omega,const int_t xmin,const int_t xmax,const int_t ymin,const int_t ymax,const int_t zmin,const int_t zmax) = 0;
    };

    template <class T>
    void Relaxation<T>::relax(Domain<T> *** domains,T omega,int_t nsmooth,int_t level,int_t num_devices){
        
        #pragma omp parallel
        #pragma omp single nowait
        #pragma omp taskgroup
        {
            for (int_t s = 0;s<smoothing_multiplier()*nsmooth;s++){
                // Swapping domains
                //#pragma omp taskgroup
                {
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    Domain<T> & domain = *domains[gpuid][level];
                    #pragma omp task default(none) shared(domain) depend(inout:domain.uprev->at[0],domain.u->at[0])
                    if (requires_duplicate_solution()) domain.swap_u();
                }
                // Applying boundary conditions
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->east->update(*domains[gpuid][level],requires_duplicate_solution());
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->west->update(*domains[gpuid][level],requires_duplicate_solution());
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->north->update(*domains[gpuid][level],requires_duplicate_solution());
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->south->update(*domains[gpuid][level],requires_duplicate_solution());
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->top->update(*domains[gpuid][level],requires_duplicate_solution());
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    domains[gpuid][level]->bottom->update(*domains[gpuid][level],requires_duplicate_solution());
                }

                // Computing halo values
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,EAST);
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,WEST);
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,NORTH);
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,SOUTH);
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,TOP);
                }
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    this->fill_send_buffer(*domains[gpuid][level],omega,BOTTOM);
                }
                // Looping over interior points
                for (int_t gpuid=0;gpuid<num_devices;gpuid++){
                    Domain<T> & domain = *domains[gpuid][level];
                    DeviceArray<T>& u = *domain.u;

                    const int_t xmin = domain.west->is_non_eliminated();
                    const int_t xmax = u.shape[0]-domain.east->is_non_eliminated();
                    const int_t ymin = domain.south->is_non_eliminated();
                    const int_t ymax = u.shape[1]-domain.north->is_non_eliminated();
                    const int_t zmin = domain.bottom->is_non_eliminated();
                    const int_t zmax = u.shape[2]-domain.top->is_non_eliminated();
                    #pragma omp task default(none) shared(domain,omega)\
                    firstprivate(xmin,xmax,ymin,ymax,zmin,zmax)\
                    depend(inout:domain.u->at[0],domain.uprev->at[0]) \
                    depend(in:domain.east->arr.at[0],domain.west->arr.at[0],domain.north->arr.at[0],\
                    domain.south->arr.at[0],domain.top->arr.at[0],domain.bottom->arr.at[0])
                    relaxation_kernel(domain,omega,xmin,xmax,ymin,ymax,zmin,zmax);
                }
                }
            }
        }
    }

    template <class T>
    void Relaxation<T>::fill_send_buffer(Domain<T>& domain, const T omega, Location_t loc)
    {
        Boundary<T> * bound;
        switch(loc){
            case NORTH:
                bound = domain.north;
                break;
            case SOUTH:
                bound = domain.south;
                break;
            case EAST:
                bound = domain.east;
                break;
            case WEST:
                bound = domain.west;
                break;
            case TOP:
                bound = domain.top;
                break;
            case BOTTOM:
                bound = domain.bottom;
                break;
        }
        if (bound->is_internal_boundary()){
            OMPBoundary<T> * omp_bound = (OMPBoundary<T> *) bound;
            #pragma omp task default(none) shared(domain,omega) firstprivate(omp_bound) depend(in:bound->arr.at[0]) depend(out:omp_bound->send_buffer.at[0])
            {
                if (requires_duplicate_solution()) {
                    omp_bound->fill_send_buffer(*domain.uprev,*domain.f,domain.settings,domain,omega);  
                }
                else {
                    omp_bound->fill_send_buffer(*domain.u,*domain.f,domain.settings,domain,omega);
                }      
            }
        }
    }
}

#endif
