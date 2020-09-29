
#include "bbn_headers.h"

double bbn_bhf_poly_smoother
(const double r,const double rmax,const double rmin)
;
double bbn_bhf_poly_smoother
(const double r,const double rmax,const double rmin)
{
const double x = (2*r-rmax-rmin)/(rmax-rmin);
return
(429.0/4096.0)*pow(x, 14) + (231.0/2048.0)*pow(x, 13) - 3003.0/4096.0*
pow(x, 12) - 819.0/1024.0*pow(x, 11) + (9009.0/4096.0)*pow(x, 10) + 
(5005.0/2048.0)*pow(x, 9) - 15015.0/4096.0*pow(x, 8) - 2145.0/512.0*
pow(x, 7) + (15015.0/4096.0)*pow(x, 6) + (9009.0/2048.0)*pow(x, 5) - 
9009.0/4096.0*pow(x, 4) - 3003.0/1024.0*pow(x, 3) + (3003.0/4096.0)*
pow(x, 2) + (3003.0/2048.0)*x + 1619.0/4096.0;
}



