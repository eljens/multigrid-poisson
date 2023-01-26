#ifndef POISSON_LIBRARY
#define POISSON_LIBRARY
#include "domain.h"
#include "injection.h"
#include "fullweighting.h"
#include "trilinearinterpolation.h"
#include "grid.h"
#include "jacobi.h"
#include "vcycle.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

using std::setw;
using std::cout;
using std::string;

namespace Poisson{
    template <class T,template<class> class R>
    class PoissonSolver {
        private:
            Domain<T> ** domains;
            Settings settings;
            Grid grid;
            R<T> restriction_type;
            TrilinearInterpolation<T> trilinearinterpolation;
            bool is_dirichlet;
            bool is_verbose = false;
            uint_t iter = 0;
            T rel_res = 0.0;
            T wtime = 0.0;
            void alloc();
        public:
            PoissonSolver(Settings & _settings,bool _is_dirichlet);
            PoissonSolver(Settings & _settings);
            ~PoissonSolver();
            void init();
            void init_zero();
            void to_device();
            void to_host();
            void verbose(bool onoff);
            void solve(const string solver="Vcycle",string smoother="Jacobi",const int_t nsmooth=10,T omega=-1.0);
            void save(string file_name="u.vtk");
            void save_all(string ufile_name="u.vtk",string ffile_name="f.vtk",string rfile_name="r.vtk");
            T relative_residual();
            T solve_time();
            T solve_iterations();
            DeviceArray<T> & get_rhs();
            const DeviceArray<T> & get_solution();
            const DeviceArray<T> & get_residual();
            Boundary<T> & get_top_bc();
            Boundary<T> & get_bottom_bc();
    };

    template <class T,template<class> class R>
    PoissonSolver<T,R>::PoissonSolver(Settings & _settings, bool _is_dirichlet) : 
        settings(_settings), 
        grid(_settings,_settings.levels), is_dirichlet(_is_dirichlet)
    {
        this->alloc();
    }

    template <class T,template<class> class R>
    PoissonSolver<T,R>::PoissonSolver(Settings & _settings) : 
        settings(_settings),
        grid(_settings,_settings.levels), is_dirichlet(false)
    {
        this->alloc();
    }

    template <class T,template<class> class R>
    PoissonSolver<T,R>::~PoissonSolver()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            delete domains[l];
        }
        delete[] domains;
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::alloc()
    {
        domains = new Domain<T>*[settings.levels];
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l] = new Domain<T>(grid.domainsettings[l],is_dirichlet);
        }
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::init()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            if (l==0){
                domains[l]->init(&ufun,&ffun,&dudxfun,&dudyfun);
            }
            else {
                domains[l]->init_zero();
            }
        }
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::init_zero()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l]->init_zero();
        }
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::to_device()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l]->to_device();
        }
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::to_host()
    {
        domains[0]->to_host();
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::verbose(bool onoff)
    {
        this->is_verbose = onoff;
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::solve(const string solver,string smoother,const int_t nsmooth,T omega){
        wtime = omp_get_wtime();
        bool use_vcycle = ((solver.compare("vcycle")==0) || (solver.compare("Vcycle")==0));
        bool use_jacobi = ((solver.compare("jacobi")==0) || (solver.compare("Jacobi")==0));
        bool use_gaussseidel = ((solver.compare("gaussseidel")==0) || (solver.compare("GaussSeidel")==0));
        if (!(use_vcycle || use_jacobi)){
            if ((solver.compare("gaussseidel")!=0) && (solver.compare("GaussSeidel")!=0)){
                throw std::invalid_argument("solver must be \"Jacobi\", \"GaussSeidel\", or \"Vcycle\" but was \""+solver+"\"");
            }
        }

        // Changing Omega dependent on method
        if (use_jacobi || (smoother.compare("jacobi")==0) || (smoother.compare("Jacobi")==0)){
            if (omega < 0) omega = 6.0/7.0;
            smoother = "jacobi";
        }
        if (use_gaussseidel || (smoother.compare("gaussseidel")==0) || (smoother.compare("GaussSeidel")==0)){
            if (omega < 0) omega = 1.05;
            smoother = "gauss_seidel";
        }

        iter = 0;
        T fnorm = domains[0]->f->infinity_norm();
        residual<T>(*domains[0]);
        rel_res = domains[0]->r->infinity_norm() / fnorm;

        if (is_verbose){
            cout << setw(4) << 0 << ": Initial residual: " << setw(8) << rel_res << endl;
        }
        for(iter = 0;iter<settings.maxiter;iter++){
            if (use_vcycle){
                Vcycle<T>(this->domains,this->restriction_type,this->trilinearinterpolation,omega,0,settings.levels,smoother,nsmooth);
            }
            else if (use_jacobi) {
                jacobi<T>(*domains[0],omega);
            }
            else {
                gaussseidel<T>(*domains[0],omega);
            }
            residual<T>(*domains[0]);
            T norm = domains[0]->r->infinity_norm() / fnorm;
            if (is_verbose){
                cout << setw(4) << iter+1 << ": Relative residual: " << setw(8) << norm << endl;
            }
            if (norm > rel_res){
                rel_res = norm;
                break;
            }
            else if (std::abs(norm-rel_res) < this->settings.tolerance){
                rel_res = norm;
                break;
            }
            rel_res = norm;
        }
        wtime = omp_get_wtime()-wtime;
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::save(string file_name)
    {
        domains[0]->save(file_name);
    }

    template <class T,template<class> class R>
    void PoissonSolver<T,R>::save_all(string ufile_name,string ffile_name,string rfile_name)
    {
        domains[0]->save_all(ufile_name,ffile_name,rfile_name);
    }

    template <class T,template<class> class R>
    T PoissonSolver<T,R>::relative_residual()
    {
        return rel_res;
    }

    template <class T,template<class> class R>
    T PoissonSolver<T,R>::solve_time()
    {
        return wtime;
    }

    template <class T,template<class> class R>
    T PoissonSolver<T,R>::solve_iterations()
    {
        return iter;
    }

    template <class T,template<class> class R>
    DeviceArray<T> & PoissonSolver<T,R>::get_rhs(){
        return *(domains[0]->f);
    }

   template <class T,template<class> class R>
    const DeviceArray<T> & PoissonSolver<T,R>::get_solution(){
        return *(domains[0]->u);
    }

    template <class T,template<class> class R>
    const DeviceArray<T> & PoissonSolver<T,R>::get_residual(){
        return *(domains[0]->r);
    }

    template <class T,template<class> class R>
    Boundary<T> & PoissonSolver<T,R>::get_top_bc(){
        return *(domains[0]->top);
    }

    template <class T,template<class> class R>
    Boundary<T> & PoissonSolver<T,R>::get_bottom_bc(){
        return *(domains[0]->bottom);
    }
}
#endif
