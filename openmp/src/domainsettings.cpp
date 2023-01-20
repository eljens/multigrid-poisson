#include "domainsettings.h"
#include <string>
#include <stdexcept>
using std::to_string;
using std::string;

DomainSettings::DomainSettings() : Settings() {
    // Nothing so far
}

DomainSettings::DomainSettings(Settings & _settings,const uint_t l) : Settings(_settings) {
    if (l>0){
        uint_t factor = 1;
        uint_t original_dims[3] = {this->dims[0],this->dims[1],this->dims[2]};
        this->dims[0] = (this->dims[0] >> l)+1;
        this->dims[1] = (this->dims[1] >> l)+1;
        this->dims[2] = (this->dims[2] >> l)+1;
        this->h *= (double_t) (factor<<l);
        if ((((this->dims[0]-1) << l)+1 != original_dims[0])
        || (((this->dims[1]-1) << l)+1 != original_dims[1])
        || (((this->dims[2]-1) << l)+1 != original_dims[2]))
        {
            string errormsg = "Invalid grid of size: ("+to_string(this->dims[0])+",";
            errormsg += to_string(this->dims[1])+","+to_string(this->dims[2])+")";
            throw std::invalid_argument(errormsg);
        }
    }
}

DomainSettings::~DomainSettings() {
    // Nothing so far
}

ostream& operator<<(ostream& os, const DomainSettings& settings){
    os << "DomainSettings:" <<endl;
    os << "\tDomain size: (" << settings.dims[0] << "," << settings.dims[1] << "," << settings.dims[2] << ")" << endl;
    os << "\tGridSpacing: " << settings.h << endl;
    return os;
}