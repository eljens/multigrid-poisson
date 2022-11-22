#include "../include/settings.h"

Settings::Settings(){
      // Does nothing
}

Settings::Settings(Settings & settings){
      this->dims[0] = settings.dims[0];
      this->dims[1] = settings.dims[1];
      this->dims[2] = settings.dims[2];
      this->maxiter = settings.maxiter;
      this->tolerance = settings.tolerance;
      this->host = settings.host;
      this->dev = settings.dev;
      this->h = settings.h;
      this->origin[0] = settings.origin[0];
      this->origin[1] = settings.origin[1];
      this->origin[2] = settings.origin[2];
}

Settings::~Settings(){

}

std::ostream& operator<<(std::ostream& os, const Settings& settings)
{
      os << "Settings:" <<endl;
      os << "\tDomain size: (" << settings.dims[0] << "," << settings.dims[1] << "," << settings.dims[2] << ")" << endl;
      os << "\tMax iterations: " << settings.maxiter << endl;
      os << "\tTolerance " << settings.tolerance << endl;
      os << "\tHost id: " << settings.host << endl;
      os << "\tDevice id: " << settings.host << endl;
      return os;
}