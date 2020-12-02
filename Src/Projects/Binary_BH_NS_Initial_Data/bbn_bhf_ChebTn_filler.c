
#include "bbn_headers.h"

void bbn_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
;
void bbn_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
{
assert(N==10);
a[0] = 
(1263.0/131072.0)*ddfddr*pow(rfill, 2) - 5359.0/65536.0*dfdr*rfill + 
(1.0/2.0)*fr0 + (1.0/2.0)*fr1;
a[1] = 
(521.0/131072.0)*ddfddr*pow(rfill, 2) - 3461.0/65536.0*dfdr*rfill - 
19845.0/32768.0*fr0 + (19845.0/32768.0)*fr1;
a[2] = 
(1.0/16384.0)*rfill*(-223*ddfddr*rfill + 1470*dfdr);
a[3] = 
-441.0/65536.0*ddfddr*pow(rfill, 2) + (2205.0/32768.0)*dfdr*rfill + 
(2205.0/16384.0)*fr0 - 2205.0/16384.0*fr1;
a[4] = 
(147.0/32768.0)*rfill*(ddfddr*rfill - 2*dfdr);
a[5] = 
(219.0/65536.0)*ddfddr*pow(rfill, 2) - 567.0/32768.0*dfdr*rfill - 
567.0/16384.0*fr0 + (567.0/16384.0)*fr1;
a[6] = 
(9.0/16384.0)*rfill*(-ddfddr*rfill + 2*dfdr);
a[7] = 
-169.0/262144.0*ddfddr*pow(rfill, 2) + (405.0/131072.0)*dfdr*rfill + 
(405.0/65536.0)*fr0 - 405.0/65536.0*fr1;
a[8] = 
(5.0/131072.0)*rfill*(ddfddr*rfill - 2*dfdr);
a[9] = 
(15.0/262144.0)*ddfddr*pow(rfill, 2) - 35.0/131072.0*dfdr*rfill - 35.0/
65536.0*fr0 + (35.0/65536.0)*fr1;
}

