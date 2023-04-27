#include "parser.h"
#include "domainsettings.h"
#include "libpoisson.h"
#include "definitions.h"
#include "devicearray.h"
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using Poisson::PoissonSolver;
using Poisson::Settings;
using Poisson::parser;
using Poisson::FullWeighting;
using Poisson::Injection;
using Poisson::GaussSeidel;
using Poisson::TrilinearInterpolation;
using Poisson::Jacobi;
using Poisson::DeviceArray;
using Poisson::int_t;
using Poisson::uint_t;

int main(int argc, char * argv[]){
    bool is_dirichlet = false;
    uint_t num_devices = omp_get_num_devices();
    Settings settings = parser(argc,argv);
    PoissonSolver<double_t,Injection,TrilinearInterpolation,Jacobi> solver(num_devices,settings,is_dirichlet);
    solver.init();
    solver.verbose(true);
    solver.to_device();

    solver.solve("vcycle",4);
    cout << "It took " << solver.solve_time() << " seconds to run ";
    cout << solver.solve_iterations() << " Fcycles"<<endl;

    solver.to_host();
    
    for (int_t gpuid = 0; gpuid < num_devices; gpuid++){
        if (settings.print_result){
            solver.save_all(gpuid,
                "results/u_"+std::to_string(gpuid)+".vtk",
                "results/f_"+std::to_string(gpuid)+".vtk",
                "results/r_"+std::to_string(gpuid)+".vtk");
        }

        const DeviceArray<double_t> & u = solver.get_solution(gpuid);

        double_t err = 0.0;

        #pragma omp parallel for collapse(3) reduction(max:err)
        for (int_t i = 0;i<u.shape[0];i++){
            for (int_t j = 0;j<u.shape[1];j++){
                for (int_t k = 0;k<u.shape[2];k++){
                    double_t u_true = 
                        Poisson::ufun(settings.origin[0]+settings.lengthx*gpuid+((double_t) i)*settings.h,
                        settings.origin[1]+((double_t) j)*settings.h,
                        settings.origin[2]+((double_t) k)*settings.h);
                    double tmp = abs(u.at[u.idx(i,j,k)]-u_true);
                    err = std::max(err,tmp);
                }
            }
        }

        cout << "Maximal error: " << setw(8) << err << endl;
    }

    return EXIT_SUCCESS;
}