
#include "bh_header.h"

void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
;
void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
{
assert(N==6);
a[0] = 
(3.0/128.0)*ddfddr*pow(rfill, 2) - 37.0/256.0*dfdr*rfill + (49.0/
128.0)*fr0 + (79.0/128.0)*fr1;
a[1] = 
-1.0/256.0*ddfddr*pow(rfill, 2) - 5.0/256.0*dfdr*rfill - 35.0/64.0*
fr0 + (35.0/64.0)*fr1;
a[2] = 
-1.0/32.0*ddfddr*pow(rfill, 2) + (11.0/64.0)*dfdr*rfill + (5.0/32.0)*
fr0 - 5.0/32.0*fr1;
a[3] = 
(3.0/512.0)*ddfddr*pow(rfill, 2) + (7.0/512.0)*dfdr*rfill + (5.0/
128.0)*fr0 - 5.0/128.0*fr1;
a[4] = 
(1.0/128.0)*ddfddr*pow(rfill, 2) - 7.0/256.0*dfdr*rfill - 5.0/128.0*
fr0 + (5.0/128.0)*fr1;
a[5] = 
-1.0/512.0*ddfddr*pow(rfill, 2) + (3.0/512.0)*dfdr*rfill + (1.0/128.0)*
fr0 - 1.0/128.0*fr1;
}

