#include "problem_definition.h"
namespace Poisson {
    #define kx 3
    #define ky 2
    #define kz 3.56

    double_t ufun(double_t x,double_t y,double_t z){
        return sin(kx*x)*sin(ky*y)*sin(kz*z);
    }

    double_t ffun(double_t x,double_t y,double_t z){
        return -(kx*kx+ky*ky+kz*kz)*sin(kx*x)*sin(ky*y)*sin(kz*z);
    }

    double_t dudxfun(double_t x,double_t y,double_t z){
        return kx*cos(kx*x)*sin(ky*y)*sin(kz*z);
    }

    double_t dudyfun(double_t x,double_t y,double_t z){
        return ky*sin(kx*x)*cos(ky*y)*sin(kz*z);
    }
}