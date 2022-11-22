#include "../include/domainsettings.h"

DomainSettings::DomainSettings() : Settings() {
    // Nothing so far
}

DomainSettings::DomainSettings(Settings & _settings,const uint_t l) : Settings(_settings) {
    if (l>0){
        uint_t factor = 2;
        for(uint_t i=0;i<l-1;i++){
            factor*=2;
        }
        this->dims[0] = ((uint_t) this->dims[0]/factor)+1;
        this->dims[1] = ((uint_t) this->dims[1]/factor)+1;
        this->dims[2] = ((uint_t) this->dims[2]/factor)+1;
        this->h *= factor;
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