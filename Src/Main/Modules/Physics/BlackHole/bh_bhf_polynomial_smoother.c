
#include "bh_header.h"

double bh_bhf_poly_smoother
(const double r,const double rmax,const double rmin)
;
double bh_bhf_poly_smoother
(const double r,const double rmax,const double rmin)
{
const double x = (2*r-rmax-rmin)/(rmax-rmin);
return
(35.0/256.0)*pow(x, 9) - 45.0/64.0*pow(x, 7) + (189.0/128.0)*
pow(x, 5) - 105.0/64.0*pow(x, 3) + (315.0/256.0)*x + 1.0/
2.0;
}

