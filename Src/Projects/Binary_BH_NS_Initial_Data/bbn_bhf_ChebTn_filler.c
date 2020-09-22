
#include "bbn_headers.h"

void bbn_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill)
;
void bbn_bhf_ChebTn_extrapolate
(double *const a,const double fr0,const double fr1,const double dfdr,const double ddfddr,const double rfill)
{
a[0] = 
(143.0/131072.0)*ddfddr*pow(rfill, 2) - 1859.0/65536.0*dfdr*rfill + 
(20899.0/32768.0)*fr0 + (11869.0/32768.0)*fr1;
a[1] = 
(143.0/131072.0)*ddfddr*pow(rfill, 2) - 2145.0/65536.0*dfdr*rfill - 
9009.0/16384.0*fr0 + (9009.0/16384.0)*fr1;
a[2] = 
-13.0/16384.0*ddfddr*pow(rfill, 2) + (91.0/8192.0)*dfdr*rfill - 819.0/
4096.0*fr0 + (819.0/4096.0)*fr1;
a[3] = 
-91.0/65536.0*ddfddr*pow(rfill, 2) + (1001.0/32768.0)*dfdr*rfill + 
(273.0/8192.0)*fr0 - 273.0/8192.0*fr1;
a[4] = 
-21.0/32768.0*ddfddr*pow(rfill, 2) + (329.0/16384.0)*dfdr*rfill + 
(567.0/8192.0)*fr0 - 567.0/8192.0*fr1;
a[5] = 
(9.0/65536.0)*ddfddr*pow(rfill, 2) + (133.0/32768.0)*dfdr*rfill + 
(189.0/8192.0)*fr0 - 189.0/8192.0*fr1;
a[6] = 
(5.0/16384.0)*ddfddr*pow(rfill, 2) - 19.0/8192.0*dfdr*rfill - 21.0/
4096.0*fr0 + (21.0/4096.0)*fr1;
a[7] = 
(41.0/262144.0)*ddfddr*pow(rfill, 2) - 239.0/131072.0*dfdr*rfill - 
207.0/32768.0*fr0 + (207.0/32768.0)*fr1;
a[8] = 
(5.0/131072.0)*ddfddr*pow(rfill, 2) - 33.0/65536.0*dfdr*rfill - 63.0/
32768.0*fr0 + (63.0/32768.0)*fr1;
a[9] = 
(1.0/262144.0)*ddfddr*pow(rfill, 2) - 7.0/131072.0*dfdr*rfill - 7.0/
32768.0*fr0 + (7.0/32768.0)*fr1;
}

