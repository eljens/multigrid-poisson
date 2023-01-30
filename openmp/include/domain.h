#ifndef TOP3D_DOMAIN
#define TOP3D_DOMAIN

#include "definitions.h"
#include "settings.h"
#include "boundary.h"
#include "neumann.h"
#include "dirichlet.h"
#include "problem_definition.h"
#include "residual.h"
#include "halo.h"
#include "restriction.h"
#include "prolongation.h"
#include "vcycle.h"

using std::cout;
using std::swap;

namespace Poisson{
	template <typename T>
	class Domain;
	
	template <typename T>
	class Relaxation;

	template <typename T>
	class Jacobi;

	template <typename T>
	class GaussSeidel;

	template <typename T>
	void residual(Domain<T>& domain);

	template <class T>
	class Domain
	{
		private:
			Halo halo;

			// Granting access to Jacobi
			friend class Relaxation<T>;

			friend class Jacobi<T>;

			friend class GaussSeidel<T>;

			friend void residual<T>(Domain<T>& domain);

			void init_f(funptr ffun);

			bool is_initialized;

		protected:
			void write_bc_to(DeviceArray<T> & uarr);
			
		public:
			Settings settings;
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

			Domain() : is_initialized(false), settings() {

			}

			Domain(Settings & _settings,bool is_dirichlet);

			~Domain();

			void init(funptr ufun,funptr ffun,funptr dudxfun,funptr dudyfun);

			void init_zero();

			void to_device();

			void to_host();

			void save(string ufilename="u.vtk");

			void save_all(string ufilename="u.vtk",string ffilename="f.vtk",string rfilename="u.vtk");

			void save_halo();

			void save_restrict_prolong(Restriction<T> & restriction,Prolongation<T> & prolongation);

			void swap_u();
	};

	template<class T>
	Domain<T>::Domain(Settings & _settings,bool is_dirichlet) : 
		halo(!is_dirichlet,!is_dirichlet,!is_dirichlet,!is_dirichlet,0,0),is_initialized(true),settings(_settings)
	{
		//cout << "Created domain with settings " << endl;
		//cout << settings;
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
		if (is_initialized){
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
		this->uprev->init_zero();
		this->uprev->init_zero();
		this->r->init_zero();
	}

	template<class T>
	void Domain<T>::init_zero(){
		this->north->init_zero();
		this->south->init_zero();
		this->east->init_zero();
		this->west->init_zero();
		this->top->init_zero();
		this->bottom->init_zero();
		this->f->init_zero();
		this->uprev->init_zero();
		this->uprev->init_zero();
		this->r->init_zero();
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
	void Domain<T>::save(string ufilename){
		this->u->print(settings,ufilename.c_str());
		cout << "Saved " << ufilename << endl;
	}

	template<class T>
	void Domain<T>::save_all(string ufilename,string ffilename,string rfilename){
		this->u->print(settings,ufilename.c_str());
		cout << "Saved " << ufilename << endl;
		this->f->print(settings,ffilename.c_str());
		cout << "Saved " << ffilename << endl;
		this->r->print(settings,rfilename.c_str());
		cout << "Saved " << rfilename << endl;
	}

	template<class T>
	void Domain<T>::save_halo(){
		//this->f->print_halo(settings,"results/f_halo.vtk");
		this->u->print_halo(settings,"results/u_halo.vtk");
		this->r->print_halo(settings,"results/r_halo.vtk");
	}

	template<class T>
	void Domain<T>::save_restrict_prolong(Restriction<T> & restriction,Prolongation<T> & prolongation){
		DeviceArray<T> u_small(settings.dev,settings.dims[0]/2+1,settings.dims[1]/2+1,settings.dims[2]/2+1);
		DeviceArray<T> u_large(settings.dev,settings.dims[0],settings.dims[1],settings.dims[2]);
		u_small.to_device();
		u_large.to_device();
		restriction.restrict_to(*(this->u),u_small);
		prolongation.prolong(u_small,u_large);
		u_small.to_host();
		u_large.to_host();
		u_small.print(settings,"results/u_restricted.vtk");
		u_large.print(settings,"results/u_prolonged.vtk");
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
		#pragma omp parallel for collapse(3) schedule(static,CHUNK_SIZE)
		for(int_t i = 0;i<this->f->shape[0];i++){
			for(int_t j = 0;j<this->f->shape[1];j++){
				for(int_t k = 0;k<this->f->shape[2];k++){
					this->f->at[this->f->idx(i,j,k)] = ffun(x0+((T)i)*h,y0+((T)j)*h,z0+((T)k)*h);
				}
			}
		}
	}
}
#endif
