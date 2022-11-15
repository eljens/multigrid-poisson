#ifndef TOP3D_DOMAIN
#define TOP3D_DOMAIN

#include "definitions.h"
#include "parser.h"
#include "boundary.h"
#include "neumann.h"
#include "dirichlet.h"
#include "problem_definition.h"
#include "jacobi.h"

using std::cout;
using std::swap;

template <typename T>
class Domain;

template <typename T>
void jacobi(Domain<T>& domain, T omega);

template <class T>
class Domain
{
	private:
		Settings & settings;

		// Granting access to Jacobi
		friend void jacobi<T>(Domain<T>& domain, T omega);

		void init_f(funptr ffun){
			const T x0 = settings.origin[0];
            const T y0 = settings.origin[1];
            const T z0 = settings.origin[2];
            const T h = settings.h;
			for(int_t i = 0;i<this->f->shape[0];i++){
                for(int_t j = 0;j<this->f->shape[1];j++){
                    for(int_t k = 0;k<this->f->shape[2];k++){
                        this->f->at[this->f->idx(i,j,k)] = ffun(x0+((T)i)*h,y0+((T)j)*h,z0+((T)k)*h);
                    }
                }
            }
		}

		void init_u(DeviceArray<T> * uarr){
			for(int_t i = 0;i<uarr->shape[0];i++){
                for(int_t j = 0;j<uarr->shape[1];j++){
                    for(int_t k = 0;k<uarr->shape[2];k++){
                        uarr->at[uarr->idx(i,j,k)] = (T) 0.0;
                    }
                }
            }
		}
	protected:
		bool even = true;

		void write_bc_to(DeviceArray<T> * uarr) {
			this->north->write_to(uarr,this->settings);
			this->south->write_to(uarr,this->settings);
			this->east->write_to(uarr,this->settings);
			this->west->write_to(uarr,this->settings);
			this->top->write_to(uarr,this->settings);
			this->bottom->write_to(uarr,this->settings);
		}
		
	public:
		DeviceArray<T> * u;
		DeviceArray<T> * uprev;
		DeviceArray<T> * f;
		Boundary<T> * north;
		Boundary<T> * south;
		Boundary<T> * east;
		Boundary<T> * west;
		Boundary<T> * top;
		Boundary<T> * bottom;
		uint_t halox[2] = {0,0};
		uint_t haloy[2] = {0,0};
		uint_t haloz[2] = {0,0};

		Domain(Settings & _settings) : settings(_settings)
		{
			cout << "Created domain with settings " << endl;
			cout << settings;
			this->u = new DeviceArray<T>(settings.dev,settings.dims[0],settings.dims[1],settings.dims[2]);
			this->uprev = new DeviceArray<T>(settings.dev,settings.dims[0],settings.dims[1],settings.dims[2]);
			this->f = new DeviceArray<T>(settings.dev,settings.dims[0],settings.dims[1],settings.dims[2]);

			this->north = new Dirichlet<T>(settings.dev,NORTH,settings.dims[0],1,settings.dims[2]);
			this->south = new Dirichlet<T>(settings.dev,SOUTH,settings.dims[0],1,settings.dims[2]);
			this->east = new Dirichlet<T>(settings.dev,EAST,1,settings.dims[1],settings.dims[2]);
			this->west = new Dirichlet<T>(settings.dev,WEST,1,settings.dims[1],settings.dims[2]);
			this->top = new Dirichlet<T>(settings.dev,TOP,settings.dims[0],settings.dims[1],1);
			this->bottom = new Dirichlet<T>(settings.dev,BOTTOM,settings.dims[0],settings.dims[1],1);
		}

		~Domain(){
			delete this->u;
			delete this->uprev;
			delete this->f;
			delete this->north;
			delete this->south;
			delete this->east;
			delete this->west;
			delete this->top;
			delete this->bottom;
		}

		void init(funptr ufun,funptr ffun,funptr dudxfun,funptr dudyfun){
			cout << "Initializing boundaries" << endl;
			this->north->init(ufun,dudxfun,dudyfun,this->settings);
			this->south->init(ufun,dudxfun,dudyfun,this->settings);
			this->east->init(ufun,dudxfun,dudyfun,this->settings);
			this->west->init(ufun,dudxfun,dudyfun,this->settings);
			this->top->init(ufun,dudxfun,dudyfun,this->settings);
			this->bottom->init(ufun,dudxfun,dudyfun,this->settings);
			cout << "Initialized boundaries" << endl;
			this->init_f(ffun);
			this->init_u(this->u);
			this->init_u(this->uprev);

			this->write_bc_to(this->u);
			this->write_bc_to(this->uprev);
		}

		void to_device(){
			this->u->to_device();
			this->uprev->to_device();
			this->f->to_device();
			this->north->to_device();
			this->south->to_device();
			this->east->to_device();
			this->west->to_device();
			this->top->to_device();
			this->bottom->to_device();
		}

		void to_host(){
			this->u->to_host();
			this->uprev->to_host();
		}

		void save_result(){
			this->f->to_vtk_file("f.vtk");
			this->u->to_vtk_file("u.vtk");
			this->uprev->to_vtk_file("uprev.vtk");
		}

		void swap_u(){
			// DeviceArray<T> * utmp = this->u;
			// this->u = this->uprev;
			// this->uprev = utmp;
			//cout << "Address of u is " << this->u << " and address of uprev is " << this->uprev << endl;
			//swap(u,&uprev);
			this->even = !(this->even);
		}
};

#endif

