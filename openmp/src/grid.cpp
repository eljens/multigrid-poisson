#include "grid.h"
#include <iostream>

using std::cout;
namespace Poisson {
    Grid::Grid(const Settings & _settings,const uint_t _levels) : 
        levels(_levels) {
            this->domainsettings = new DomainSettings[_levels];
            for(uint_t l=0;l<_levels;l++){
                this->domainsettings[l] = DomainSettings(_settings,l);
                cout << this->domainsettings[l];
            }
    }

    Grid::~Grid(){
        delete[] this->domainsettings;
    }
}