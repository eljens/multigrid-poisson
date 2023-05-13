#ifndef POISSON_LIBRARY
#define POISSON_LIBRARY
#include "domain.h"
#include "injection.h"
#include "fullweighting.h"
#include "trilinearinterpolation.h"
#include "grid.h"
#include "jacobi.h"
#include "gaussseidel.h"
#include "vcycle.h"
#include "fcycle.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <fstream>
#include "boundary.h"

using std::ofstream;
using std::ios;
using std::setw;
using std::cout;
using std::string;
using Poisson::DIRICHLET;
using Poisson::NEUMANN;

namespace Poisson{
    template <class T,template<class> class R,template<class> class P,template<class> class S>
    class PoissonSolver {
        private:
            uint_t num_devices;
            Domain<T> *** domains;
            Settings settings;
            BoundaryCondition BC;
            Grid grid;
            R<T> restriction_type;
            P<T> prolongation;
            S<T> relaxation;
            bool is_verbose = false;
            uint_t iter = 0;
            double_t rel_res = 0.0;
            double_t wtime = 0.0;
            void alloc();
        public:
            PoissonSolver(uint_t ndev, Settings & _settings,BoundaryCondition & _BC);
            PoissonSolver(uint_t ndev, Settings & _settings);
            ~PoissonSolver();
            void init();
            void init_zero();
            void to_device();
            void to_host();
            void verbose(bool onoff);
            void solve(const string solver="Vcycle",const int_t nsmooth=10,T omega=-1.0,string ofile="",double maxtime = 20*60);
            void save(uint_t gpuid,string file_name="u.vtk");
            void save_all(uint_t gpuid,string ufile_name="u.vtk",string ffile_name="f.vtk",string rfile_name="r.vtk");
            T relative_residual();
            T solve_time();
            T solve_iterations();
            DeviceArray<T> & get_rhs(uint_t gpuid);
            const DeviceArray<T> & get_solution(uint_t gpuid);
            const DeviceArray<T> & get_residual(uint_t gpuid);
            Boundary<T> & get_top_bc(uint_t gpuid);
            Boundary<T> & get_bottom_bc(uint_t gpuid);
    };

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    PoissonSolver<T,R,P,S>::PoissonSolver(uint_t ndev, Settings & _settings, BoundaryCondition & _BC) : 
        num_devices(ndev), settings(_settings), BC(_BC), grid(_settings,_settings.levels)
    {
        this->alloc();
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    PoissonSolver<T,R,P,S>::PoissonSolver(uint_t ndev, Settings & _settings) : 
        num_devices(ndev), settings(_settings),
        BC({NEUMANN,NEUMANN,NEUMANN,NEUMANN,DIRICHLET,DIRICHLET}), grid(_settings,_settings.levels)
    {
        this->alloc();
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    PoissonSolver<T,R,P,S>::~PoissonSolver()
    {
        for (uint_t gpuid = 0;gpuid<num_devices;gpuid++){
            for (uint_t l = 0;l<settings.levels;l++){
                delete domains[gpuid][l];
            }
            delete domains[gpuid];
        }
        delete[] domains;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::alloc()
    {
        domains = new Domain<T>**[num_devices];
        for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
            domains[gpuid] = new Domain<T>*[settings.levels];
            for (uint_t l = 0;l<settings.levels;l++){
                grid.domainsettings[l].dev = gpuid;
                domains[gpuid][l] = new Domain<T>(grid.domainsettings[l],BC,relaxation.requires_duplicate_solution(),gpuid,num_devices);
                /*// Naive slab decomposition in x dimension 
                grid.domainsettings[l].origin[0] += grid.domainsettings[l].lengthx;
                if (gpuid > 0){
                    domains[gpuid][l]->west->link(domains[gpuid-1][l]->east);
                    domains[gpuid-1][l]->east->link(domains[gpuid][l]->west);
                }*/
                //Naive slab decomposition in y dimension 
                grid.domainsettings[l].origin[1] += (grid.domainsettings[l].dims[1]-1)*grid.domainsettings[l].h;
                if (gpuid > 0){
                    domains[gpuid][l]->south->link(domains[gpuid-1][l]->north);
                    domains[gpuid-1][l]->north->link(domains[gpuid][l]->south);
                }
                /*// Naive slab decomposition in z dimension 
                grid.domainsettings[l].origin[2] += (grid.domainsettings[l].dims[2]-1)*grid.domainsettings[l].h;
                if (gpuid > 0){
                    domains[gpuid][l]->bottom->link(domains[gpuid-1][l]->top);
                    domains[gpuid-1][l]->top->link(domains[gpuid][l]->bottom);
                }*/
            }
        }
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::init()
    {   
        for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
            for (uint_t l = 0;l<settings.levels;l++){
                if (l==0){
                    domains[gpuid][l]->init(&ufun,&ffun,&dudxfun,&dudyfun);
                }
                else {
                    domains[gpuid][l]->init_zero();
                }
            }
        }
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::init_zero()
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
                    for (uint_t l = 0;l<settings.levels;l++){
                        #pragma omp task
                        {
                            domains[gpuid][l]->init_zero();
                        }
                    }
                }
            }
        }
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::to_device()
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
                    for (uint_t l = 0;l<settings.levels;l++){
                        #pragma omp task
                        {
                            domains[gpuid][l]->to_device();
                        }
                    }
                }
            }
        }
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::to_host()
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
                    for (uint_t l = 0;l<settings.levels;l++){
                        #pragma omp task
                        {
                            domains[gpuid][l]->to_host();
                        }
                    }
                }
            }
        }
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::verbose(bool onoff)
    {
        this->is_verbose = onoff;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::solve(const string solver, const int_t nsmooth,T omega,string ofile,double maxtime){
        wtime = omp_get_wtime();
        bool use_vcycle = ((solver.compare("vcycle")==0) || (solver.compare("Vcycle")==0));
        bool use_fcycle = ((solver.compare("fcycle")==0) || (solver.compare("Fcycle")==0));
        bool use_relaxation = ((solver.compare("relaxation")==0) || (solver.compare("Relaxation")==0));
        if (!(use_vcycle || use_relaxation || use_fcycle)){
            throw std::invalid_argument("solver must be \"Relaxation\", \"Vcycle\" of \"Fcycle\" but was \""+solver+"\"");
        }

        // Changing Omega dependent on method
        if (omega < 0) omega = relaxation.default_omega();

        iter = 0;
        double_t fnorm = 0.0;
        double_t rnorm = 0.0;

        // Calculating initial relative residual
        #pragma omp parallel reduction(+:fnorm) reduction(+:rnorm)
        {
            #pragma omp single nowait
            {
                #pragma omp taskgroup task_reduction(+:rnorm) task_reduction(+:fnorm)
                {
                    for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){

                        #pragma omp task default(none) shared(domains) firstprivate(gpuid) depend(inout:domains[gpuid][0]->f) in_reduction(+:fnorm)
                        {
                            double_t fnorm_task = domains[gpuid][0]->f->infinity_norm();
                            #pragma omp atomic
                            fnorm += fnorm_task;
                        }
                        #pragma omp task depend(out:domains[gpuid][0]->r)
                        {
                            residual<T>(*domains[gpuid][0]);
                        }

                        #pragma omp task default(none) shared(domains) firstprivate(gpuid) depend(in:domains[gpuid][0]->r) in_reduction(+:rnorm)
                        {
                            double_t rnorm_task = domains[gpuid][0]->r->infinity_norm();
                            #pragma omp atomic
                            rnorm += rnorm_task;
                        }
                    }
                }
            }
        }

        rel_res = rnorm / fnorm;

        if (is_verbose){
            cout << setw(4) << 0 << ": Initial residual: " << setw(8) << rel_res << endl;
        }
        ofstream out(ofile,ios::app);
        if (ofile.compare("")!=0){
            out << "#    seconds     rel_res" << endl;
        }
        for(iter = 0;iter<settings.maxiter;iter++){
            if (use_vcycle){
                for (uint_t gpuid = 0; gpuid < num_devices; gpuid++){
                    Vcycle<T>(this->domains,this->restriction_type,this->prolongation,this->relaxation,omega,0,settings.levels,num_devices,nsmooth);
                }
            }
            else if (use_fcycle){
                //for (uint_t gpuid = 0; gpuid < num_devices; gpuid++){
                //    Fcycle<T>(this->domains[gpuid],this->restriction_type,this->prolongation,this->relaxation,omega,0,settings.levels,nsmooth);
                //}
            }
            else {
                //for (uint_t gpuid = 0; gpuid<num_devices;gpuid++){
                //    relaxation.relax(*domains[gpuid][0],omega);
                //}
            }
            T norm = 0.0;
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    #pragma omp taskgroup task_reduction(+:norm)
                    {
                        for (uint_t gpuid = 0; gpuid < num_devices; gpuid++){
                            #pragma omp task default(none) shared(domains) firstprivate(gpuid) depend(out:domains[gpuid][0]->r)
                            {
                                residual<T>(*domains[gpuid][0]);
                            }
                            #pragma omp task default(none) shared(domains) firstprivate(gpuid) depend(in:domains[gpuid][0]->r) in_reduction(+:norm)
                            {   
                                double_t norm_task = domains[gpuid][0]->r->infinity_norm();
                                #pragma omp atomic
                                norm += norm_task;
                            }
                        }
                    }
                }
            }
            norm /= fnorm;
            if (ofile.compare("")!=0){
                out << setw(12) << omp_get_wtime()-wtime;
                out << setw(12) << this->relative_residual();
                out << endl;
            }
            if (is_verbose){
                cout << setw(4) << iter+1 << ": Relative residual: " << setw(8) << norm << endl;
            }
            if (iter >= settings.miniter){
                if (norm > 2.0*rel_res){
                    rel_res = norm;
                    break;
                }
                else if (std::abs(norm-rel_res) < this->settings.tolerance){
                    rel_res = norm;
                    break;
                }
                if ((omp_get_wtime()-wtime) > maxtime ){
                    cout << "WARNING: Solver reached maximum wall time without converging!" << endl;
                    break;
                }
            }
            rel_res = norm;
        }
        wtime = omp_get_wtime()-wtime;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::save(uint_t gpuid, string file_name)
    {
        domains[gpuid][0]->save(file_name);
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    void PoissonSolver<T,R,P,S>::save_all(uint_t gpuid, string ufile_name,string ffile_name,string rfile_name)
    {
        domains[gpuid][0]->save_all(ufile_name,ffile_name,rfile_name);
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    T PoissonSolver<T,R,P,S>::relative_residual()
    {
        return rel_res;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    T PoissonSolver<T,R,P,S>::solve_time()
    {
        return wtime;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    T PoissonSolver<T,R,P,S>::solve_iterations()
    {
        return iter;
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    DeviceArray<T> & PoissonSolver<T,R,P,S>::get_rhs(uint_t gpuid){
        return *(domains[gpuid][0]->f);
    }

   template <class T,template<class> class R,template<class> class P,template<class> class S>
    const DeviceArray<T> & PoissonSolver<T,R,P,S>::get_solution(uint_t gpuid){
        return *(domains[gpuid][0]->u);
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    const DeviceArray<T> & PoissonSolver<T,R,P,S>::get_residual(uint_t gpuid){
        return *(domains[gpuid][0]->r);
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    Boundary<T> & PoissonSolver<T,R,P,S>::get_top_bc(uint_t gpuid){
        return *(domains[gpuid][0]->top);
    }

    template <class T,template<class> class R,template<class> class P,template<class> class S>
    Boundary<T> & PoissonSolver<T,R,P,S>::get_bottom_bc(uint_t gpuid){
        return *(domains[gpuid][0]->bottom);
    }
}
#endif
