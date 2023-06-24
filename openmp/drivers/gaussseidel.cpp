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
using Poisson::BoundaryCondition;

int main(int argc, char * argv[]){
    uint_t num_devices = omp_get_num_devices();
    Settings settings = parser(argc,argv);

    // Defining problem
    BoundaryCondition BC;
    BC.east = NEUMANN;
    BC.west = NEUMANN;
    BC.north = NEUMANN;
    BC.south = NEUMANN;
    BC.top = DIRICHLET;
    BC.bottom = DIRICHLET;
    int_t dev_shape[3] = {(int_t) num_devices,1,1};
    if (num_devices == 8){
        dev_shape[0] = 2;
        dev_shape[1] = 2;
        dev_shape[2] = 2;
    }
    if (num_devices == 6){
        dev_shape[0] = 3;
        dev_shape[1] = 2;
        dev_shape[2] = 1;
    }
    if (num_devices == 4){
        dev_shape[0] = 1;
        dev_shape[1] = 2;
        dev_shape[2] = 2;
    }

    cout << "MG running on " << omp_get_max_threads() << " threads and " << num_devices << " devices" << endl;

    PoissonSolver<double_t,FullWeighting,TrilinearInterpolation,GaussSeidel> solver(dev_shape,settings,BC);
    solver.init();
    solver.verbose(true);
    solver.to_device();

    solver.solve("vcycle",6,0.7);
    cout << "It took " << solver.solve_time() << " seconds to run ";
    cout << solver.solve_iterations() << " Fcycles"<<endl;

    solver.to_host();

    double_t maxerr = 0.0;
    int_t gpuid = 0;
    for (int_t ii=0;ii<dev_shape[0];ii++){
        for (int_t jj=0;jj<dev_shape[1];jj++){
            for (int_t kk=0;kk<dev_shape[2];kk++){
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
                                Poisson::ufun(
                                settings.origin[0]+((double_t) i)*settings.h+ii*(settings.dims[0]-1)*settings.h,
                                settings.origin[1]+((double_t) j)*settings.h+jj*(settings.dims[1]-1)*settings.h,
                                settings.origin[2]+((double_t) k)*settings.h+kk*(settings.dims[2]-1)*settings.h);
                            double tmp = abs(u.at[u.idx(i,j,k)]-u_true);
                            err = std::max(err,tmp);
                        }
                    }
                }
                maxerr = std::max(maxerr,err);
                gpuid++;
            }
        }
    }

    cout << "Maximal error: " << setw(8) << maxerr << endl;

    if (settings.write_final_stats){
        ofstream out(settings.stats_file, ios::app);
        out << "#    seconds     abs_err     rel_res domain_size     spacing     maxiter     lengthx      levels     threads     devices" << endl;
        out << setw(12) << solver.solve_time();
        out << setw(12) << maxerr;
        out << setw(12) << solver.relative_residual();
        out << setw(12) << settings.dims[0]*settings.dims[1]*settings.dims[2]*omp_get_num_devices();
        out << setw(12) << settings.h;
        out << setw(12) << solver.solve_iterations();
        out << setw(12) << settings.lengthx;
        out << setw(12) << settings.levels;
        out << setw(12) << omp_get_max_threads();
        out << setw(12) << omp_get_num_devices();
        out << endl;
    }


    return EXIT_SUCCESS;
}