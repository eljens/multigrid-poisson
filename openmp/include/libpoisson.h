#ifndef POISSON_LIBRARY
#define POISSON_LIBRARY
#include "domain.h"
#include "injection.h"
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
    template <class T>
    class PoissonSolver {
        private:
            const T omega = 4.5/5.0;
            Domain<T> ** domains;
            Settings settings;
            Grid grid;
            Injection<T> injection;
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
            void solve(T tol=1e-5,string solver="Vcycle");
            void save(string file_name="u.vtk");
            T relative_residual();
            T solve_time();
            T solve_iterations();
    };

    template<class T>
    PoissonSolver<T>::PoissonSolver(Settings & _settings, bool _is_dirichlet) : 
        settings(_settings), 
        grid(_settings,_settings.levels), is_dirichlet(_is_dirichlet)
    {
        this->alloc();
    }

    template<class T>
    PoissonSolver<T>::PoissonSolver(Settings & _settings) : 
        settings(_settings),
        grid(_settings,_settings.levels), is_dirichlet(false)
    {
        this->alloc();
    }

    template<class T>
    PoissonSolver<T>::~PoissonSolver()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            delete domains[l];
        }
        delete[] domains;
    }

    template<class T>
    void PoissonSolver<T>::alloc()
    {
        domains = new Domain<T>*[settings.levels];
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l] = new Domain<T>(grid.domainsettings[l],is_dirichlet);
        }
    }

    template<class T>
    void PoissonSolver<T>::init()
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

    template<class T>
    void PoissonSolver<T>::init_zero()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l]->init_zero();
        }
    }

    template<class T>
    void PoissonSolver<T>::to_device()
    {
        for (uint_t l = 0;l<settings.levels;l++){
            domains[l]->to_device();
        }
    }

    template<class T>
    void PoissonSolver<T>::to_host()
    {
        domains[0]->to_host();
    }

    template<class T>
    void PoissonSolver<T>::verbose(bool onoff)
    {
        this->is_verbose = onoff;
    }

    template<class T>
    void PoissonSolver<T>::solve(T tol,string solver){
        wtime = omp_get_wtime();
        iter = 0;
        T fnorm = domains[0]->f->infinity_norm();
        rel_res = domains[0]->r->infinity_norm() / fnorm;
        bool use_vcycle = ((solver.compare("vcycle")==0) || (solver.compare("Vcycle")==0));
        if (!use_vcycle){
            if ((solver.compare("jacobi")!=0) && (solver.compare("Jacobi")==0)){
                throw std::invalid_argument("solver must be \"Jacobi\" or \"Vcycle\"");
            }
        }
        for(iter = 0;iter<settings.maxiter;iter++){
            if (use_vcycle){
                Vcycle<T>(domains,injection,trilinearinterpolation,omega,0,settings.levels);
            }
            else {
                jacobi<T>(*domains[0],omega);
            }
            residual<T>(*domains[0]);
            T norm = domains[0]->r->infinity_norm() / fnorm;
            if (is_verbose){
                cout << setw(4) << iter+1 << ": Relative residual: " << setw(8) << norm << endl;
            }
            if (std::abs(norm-rel_res) < tol){
                rel_res = norm;
                break;
            }
            rel_res = norm;
        }
        wtime = omp_get_wtime()-wtime;
    }

    template<class T>
    void PoissonSolver<T>::save(string file_name)
    {
        domains[0]->save(file_name);
    }

    template<class T>
    T PoissonSolver<T>::relative_residual()
    {
        return rel_res;
    }

    template<class T>
    T PoissonSolver<T>::solve_time()
    {
        return wtime;
    }

    template<class T>
    T PoissonSolver<T>::solve_iterations()
    {
        return iter;
    }
}
#endif
