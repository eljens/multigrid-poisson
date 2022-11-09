#ifndef TOP3D_DOMAIN
#define TOP3D_DOMAIN

#include definitions.h

class Domain
{
	private:
		uint_t nelx;
		uint_t nely;
		uint_t nelz;
	public:
		Domain(uint_t x,uint_t y,uint_t z)
		{
			this->nelx = x;
			this->nely = y;
			this->nelz = z;
		}	
};

#endif

