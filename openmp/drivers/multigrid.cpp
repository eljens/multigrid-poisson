#include "parser.h"
#include "domain.h"
#include "definitions.h"
#include "problem_definition.h"
#include "jacobi.h"
#include "residual.h"
#include "injection.h"
#include "trilinearinterpolation.h"
#include "fullweighting.h"
#include "grid.h"
#include "fcycle.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

using std::cout;
using std::endl;
using std::ofstream;
using std::ios;
using std::setw;

using Poisson::Settings;
using Poisson::Grid;
using Poisson::Domain;
using Poisson::Injection;
using Poisson::TrilinearInterpolation;
using Poisson::Jacobi;
using Poisson::Fcycle;
using Poisson::residual;
using Poisson::parser;
using Poisson::int_t;
using Poisson::uint_t;
using Poisson::double_t;
using Poisson::ufun;
using Poisson::ffun;
using Poisson::dudxfun;
using Poisson::dudyfun;

int main(int argc, char * argv[]){
    bool is_dirichlet = false;
    const Settings settings = parser(argc,argv);
    cout << settings << endl;
    Grid grid(settings,settings.levels);
    const double_t omega = 6.0/7.0;

    // Selecting relaxation type
    Poisson::Jacobi<double_t> relax;

    // Making array of domains
    Domain<double_t> * domains[settings.levels];
    for (uint_t l = 0;l<settings.levels;l++){
        domains[l] = new Domain<double_t>(grid.domainsettings[l],is_dirichlet,relax.requires_duplicate_solution());
        if (l==0){
            domains[l]->init(&ufun,&ffun,&dudxfun,&dudyfun);
        }
        else {
            domains[l]->init_zero();
        }
        domains[l]->to_device();
    }

    double_t fnorm = domains[0]->f->infinity_norm();

    Poisson::Injection<double_t> restriction;
    TrilinearInterpolation<double_t> trilinearinterpolation;
    double_t start = omp_get_wtime();
    uint_t iter = 0;
    double_t rel_res = domains[0]->r->infinity_norm() / fnorm;
    for(iter = 0;iter<settings.maxiter;iter++){
        Fcycle<double_t>(domains,restriction,trilinearinterpolation,relax,omega,0,settings.levels,4);
        residual<double_t>(*domains[0]);
        double_t norm = domains[0]->r->infinity_norm() / fnorm;
        cout << setw(4) << iter+1 << ": Relative residual: " << setw(8) << norm << endl;
        if (std::abs(norm-rel_res) < 1e-8){
            rel_res = norm;
            break;
        }
        rel_res = norm;
    }
    double_t time_taken = omp_get_wtime()-start;
    cout << endl;
    cout << "It took " << time_taken << " seconds to run " << iter << " Vcycles"  <<endl;

    residual<double_t>(*domains[0]);

    domains[0]->to_host();
    if (settings.print_result){
        domains[0]->save("results/u.vtk");
    }

    double_t err = -1.0;

    #pragma omp parallel for collapse(3) reduction(max:err)
    for (int_t i = 0;i<domains[0]->u->shape[0];i++){
        for (int_t j = 0;j<domains[0]->u->shape[1];j++){
            for (int_t k = 0;k<domains[0]->u->shape[2];k++){
                double_t u_true = 
                    ufun(settings.origin[0]+((double_t) i)*settings.h,
                    settings.origin[1]+((double_t) j)*settings.h,
                    settings.origin[2]+((double_t) k)*settings.h);
                double tmp = abs(domains[0]->u->at[domains[0]->u->idx(i,j,k)]-u_true);
                err = std::max(err,tmp);
                domains[0]->u->at[domains[0]->u->idx(i,j,k)] = u_true;
            }
        }
    }
    if (settings.print_result){
        domains[0]->save("results/u_true.vtk");
    }

    for (uint_t l = 0;l<settings.levels;l++){
        delete domains[l];
    }

    if (settings.write_final_stats){
        ofstream out("results/"+settings.stats_file, ios::app);
        out << "#    seconds     abs_err     rel_res domain_size     spacing     maxiter     lengthx      levels" << endl;
        out << setw(12) << time_taken;
        out << setw(12) << err;
        out << setw(12) << rel_res;
        out << setw(12) << settings.dims[0]*settings.dims[1]*settings.dims[2];
        out << setw(12) << settings.h;
        out << setw(12) << iter;
        out << setw(12) << settings.lengthx;
        out << setw(12) << settings.levels;
        out << endl;
    }

    cout << "Maximal error: " << setw(8) << err << endl;

    return EXIT_SUCCESS;
}