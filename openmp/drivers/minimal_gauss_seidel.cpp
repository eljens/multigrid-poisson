#include "parser.h"
#include "domainsettings.h"
#include "libpoisson.h"
#include <iostream>

using std::cout;
using std::endl;
using Poisson::PoissonSolver;
using Poisson::Settings;
using Poisson::parser;
using Poisson::Injection;
using Poisson::TrilinearInterpolation;
using Poisson::GaussSeidel;

int main(int argc, char * argv[]){
    bool is_dirichlet = false;
    Settings settings = parser(argc,argv);
    PoissonSolver<double_t,Injection,TrilinearInterpolation,GaussSeidel> solver(settings,is_dirichlet);

    solver.init();
    solver.verbose(true);
    solver.to_device();

    solver.solve("relaxation");
    cout << "It took " << solver.solve_time() << " seconds to run ";
    cout << solver.solve_iterations() << " Gauss-Seidel iterations"  <<endl;
    
    solver.to_host();
    solver.save("results/u_gaussseidel.vtk");

    return EXIT_SUCCESS;
}