#ifndef TOP3D_DOMAIN
#define TOP3D_DOMAIN

#include "definitions.h"
#include "parser.h"
#include "boundary.h"
#include "neumann.h"
#include "dirichlet.h"
#include "problem_definition.h"
#include "jacobi.h"
#include "residual.h"
#include "halo.h"

using std::cout;
using std::swap;

template <typename T>
class Domain;

template <typename T>
void jacobi(Domain<T>& domain, T omega);

template <typename T>
void residual(Domain<T>& domain);

template <class T>
class Domain
{
	private:
		Settings & settings;
		Halo halo;

		// Granting access to Jacobi
		friend void jacobi<T>(Domain<T>& domain, T omega);

		friend void residual<T>(Domain<T>& domain);

		void init_f(funptr ffun);

		void init_u(DeviceArray<T> * uarr);
	protected:
		void write_bc_to(DeviceArray<T> & uarr);
		
	public:
		DeviceArray<T> * u;
		DeviceArray<T> * uprev;
		DeviceArray<T> * f;
		DeviceArray<T> * r;
		Boundary<T> * north;
		Boundary<T> * south;
		Boundary<T> * east;
		Boundary<T> * west;
		Boundary<T> * top;
		Boundary<T> * bottom;

		Domain(Settings & _settings,bool is_dirichlet);

		~Domain();

		void init(funptr ufun,funptr ffun,funptr dudxfun,funptr dudyfun);

		void to_device();

		void to_host();

		void save();

		void save_halo();

		void swap_u();
};

template<class T>
Domain<T>::Domain(Settings & _settings,bool is_dirichlet) : 
	settings(_settings), halo(!is_dirichlet,!is_dirichlet,!is_dirichlet,!is_dirichlet,0,0)
{
	cout << "Created domain with settings " << endl;
	cout << settings;
	this->u = new DeviceArray<T>(settings,halo);
	this->uprev = new DeviceArray<T>(settings,halo);
	this->f = new DeviceArray<T>(settings,halo);
	this->r = new DeviceArray<T>(settings,halo);

	this->top = new Dirichlet<T>(settings.dev,TOP,settings.dims[0],settings.dims[1],1);
	this->bottom = new Dirichlet<T>(settings.dev,BOTTOM,settings.dims[0],settings.dims[1],1);
	if (is_dirichlet){
		this->east = new Dirichlet<T>(settings.dev,EAST,1,settings.dims[1],settings.dims[2]);
		this->west = new Dirichlet<T>(settings.dev,WEST,1,settings.dims[1],settings.dims[2]);
		this->north = new Dirichlet<T>(settings.dev,NORTH,settings.dims[0],1,settings.dims[2]);
		this->south = new Dirichlet<T>(settings.dev,SOUTH,settings.dims[0],1,settings.dims[2]);
	}
	else {
		this->east = new Neumann<T>(settings.dev,EAST,1,settings.dims[1],settings.dims[2]);
		this->west = new Neumann<T>(settings.dev,WEST,1,settings.dims[1],settings.dims[2]);
		this->north = new Neumann<T>(settings.dev,NORTH,settings.dims[0],1,settings.dims[2]);
		this->south = new Neumann<T>(settings.dev,SOUTH,settings.dims[0],1,settings.dims[2]);
	}
}

template<class T>
Domain<T>::~Domain(){
	delete this->u;
	delete this->uprev;
	delete this->f;
	delete this->r;
	delete this->north;
	delete this->south;
	delete this->east;
	delete this->west;
	delete this->top;
	delete this->bottom;
}

template<class T>
void Domain<T>::init(funptr ufun,funptr ffun,funptr dudxfun,funptr dudyfun){
	this->north->init(ufun,dudxfun,dudyfun,this->settings);
	this->south->init(ufun,dudxfun,dudyfun,this->settings);
	this->east->init(ufun,dudxfun,dudyfun,this->settings);
	this->west->init(ufun,dudxfun,dudyfun,this->settings);
	this->top->init(ufun,dudxfun,dudyfun,this->settings);
	this->bottom->init(ufun,dudxfun,dudyfun,this->settings);
	this->init_f(ffun);
	this->init_u(this->u);
	this->init_u(this->uprev);
	this->init_u(this->r);
}

template<class T>
void Domain<T>::to_device(){
	this->u->to_device();
	this->uprev->to_device();
	this->f->to_device();
	this->r->to_device();
	this->north->to_device();
	this->south->to_device();
	this->east->to_device();
	this->west->to_device();
	this->top->to_device();
	this->bottom->to_device();

	this->write_bc_to(*(this->u));
	this->write_bc_to(*(this->uprev));
}

template<class T>
void Domain<T>::to_host(){
	this->r->to_host();
	this->u->to_host();
	this->uprev->to_host();
}

template<class T>
void Domain<T>::save(){
	this->f->print(settings,"results/f.vtk");
	this->u->print(settings,"results/u.vtk");
	this->r->print(settings,"results/r.vtk");
}

template<class T>
void Domain<T>::save_halo(){
	this->f->print_halo(settings,"results/f_halo.vtk");
	this->u->print_halo(settings,"results/u_halo.vtk");
	this->r->print_halo(settings,"results/r_halo.vtk");
}

template<class T>
void Domain<T>::swap_u(){
	swap(this->u,this->uprev);
}

template<class T>
void Domain<T>::write_bc_to(DeviceArray<T> & uarr) {
	this->north->write_to(uarr,this->settings);
	this->south->write_to(uarr,this->settings);
	this->east->write_to(uarr,this->settings);
	this->west->write_to(uarr,this->settings);
	this->top->write_to(uarr,this->settings);
	this->bottom->write_to(uarr,this->settings);
}

template<class T>
void Domain<T>::init_f(funptr ffun){
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

template<class T>
void Domain<T>::init_u(DeviceArray<T> * uarr){
	for(int_t i = 0;i<uarr->shape[0]+uarr->halo.east+uarr->halo.west;i++){
        for(int_t j = 0;j<uarr->shape[1]+uarr->halo.north+uarr->halo.south;j++){
            for(int_t k = 0;k<uarr->shape[2]+uarr->halo.top+uarr->halo.bottom;k++){
                uarr->at[uarr->idx_halo(i,j,k)] = (T) 0.0;
            }
        }
    }
}

#endif

