
#include "bh_header.h"

void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
;
void bh_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill,const Uint N)
{
assert(N==10);
a[0] = 
(311.0/16384.0)*ddfddr*pow(rfill, 2) - 7809.0/65536.0*dfdr*rfill + 
(14179.0/32768.0)*fr0 + (18589.0/32768.0)*fr1;
a[1] = 
(69.0/32768.0)*ddfddr*pow(rfill, 2) - 2971.0/65536.0*dfdr*rfill - 
4851.0/8192.0*fr0 + (4851.0/8192.0)*fr1;
a[2] = 
-117.0/4096.0*ddfddr*pow(rfill, 2) + (1225.0/8192.0)*dfdr*rfill + 
(441.0/4096.0)*fr0 - 441.0/4096.0*fr1;
a[3] = 
-49.0/16384.0*ddfddr*pow(rfill, 2) + (1715.0/32768.0)*dfdr*rfill + 
(441.0/4096.0)*fr0 - 441.0/4096.0*fr1;
a[4] = 
(49.0/4096.0)*ddfddr*pow(rfill, 2) - 637.0/16384.0*dfdr*rfill - 441.0/
8192.0*fr0 + (441.0/8192.0)*fr1;
a[5] = 
(11.0/16384.0)*ddfddr*pow(rfill, 2) - 217.0/32768.0*dfdr*rfill - 63.0/
4096.0*fr0 + (63.0/4096.0)*fr1;
a[6] = 
-11.0/4096.0*ddfddr*pow(rfill, 2) + (79.0/8192.0)*dfdr*rfill + (63.0/
4096.0)*fr0 - 63.0/4096.0*fr1;
a[7] = 
(19.0/65536.0)*ddfddr*pow(rfill, 2) - 85.0/131072.0*dfdr*rfill - 9.0/
16384.0*fr0 + (9.0/16384.0)*fr1;
a[8] = 
(5.0/16384.0)*ddfddr*pow(rfill, 2) - 75.0/65536.0*dfdr*rfill - 63.0/
32768.0*fr0 + (63.0/32768.0)*fr1;
a[9] = 
-5.0/65536.0*ddfddr*pow(rfill, 2) + (35.0/131072.0)*dfdr*rfill + (7.0/
16384.0)*fr0 - 7.0/16384.0*fr1;
}

