#ifndef TOP3D_DOMAIN
#define TOP3D_DOMAIN

#include "definitions.h"
#include "parser.h"
#include "boundary.h"
#include "problem_definition.h"

using std::cout;

template <class T>
class Domain
{
	private:
		Settings & settings;
	protected:
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
			this->u = new DeviceArray<T>(settings.dev,{settings.dims[0],settings.dims[1],settings.dims[2]});
			this->uprev = new DeviceArray<T>(settings.dev,{settings.dims[0],settings.dims[1],settings.dims[2]});
			this->f = new DeviceArray<T>(settings.dev,{settings.dims[0],settings.dims[1],settings.dims[2]});

			this->north = new Dirichlet<T>(settings.dev,NORTH,{settings.dims[0],1,settings.dims[2]});
			this->south = new Dirichlet<T>(settings.dev,SOUTH,{settings.dims[0],1,settings.dims[2]});
			this->east = new Dirichlet<T>(settings.dev,EAST,{1,settings.dims[1],settings.dims[2]});
			this->west = new Dirichlet<T>(settings.dev,WEST,{1,settings.dims[1],settings.dims[2]});
			this->top = new Dirichlet<T>(settings.dev,TOP,{settings.dims[0],settings.dims[1],1});
			this->bottom = new Dirichlet<T>(settings.dev,BOTTOM,{settings.dims[0],settings.dims[1],1});
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
			this->north->init(ufun,dudxfun,dudyfun,this->settings);
			this->south->init(ufun,dudxfun,dudyfun,this->settings);
			this->east->init(ufun,dudxfun,dudyfun,this->settings);
			this->west->init(ufun,dudxfun,dudyfun,this->settings);
			this->top->init(ufun,dudxfun,dudyfun,this->settings);
			this->bottom->init(ufun,dudxfun,dudyfun,this->settings);

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
		}

		void save_result(){
			this->u->to_vtk_file();
		}
};

#endif

