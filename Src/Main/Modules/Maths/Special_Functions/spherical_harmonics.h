#include "core_lib.h"
#include "maths_general_lib.h"
#include <complex.h>

double 
d_associated_legendre_dtheta
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  );

double 
associated_legendre
  (
  const int l/* l in P^l_m(x) */,
  const int m/* m in P^l_m(x), 0 <= m <= l*/,
  const double x/* x in P^l_m(x), -1 <= x=cos(theta) <= 1 */
  );

double complex 
Ylm
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );

double complex 
dYlm_dphi
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );

double complex 
dYlm_dtheta
  (
  const int l/* l in Y_l^m(theta,phi) */, 
  const int m/* m in Y_l^m(theta,phi), -l <= m <= l */, 
  const double theta/* theta in Y_l^m(theta,phi) */,
  const double phi/*  phi in Y_l^m(theta,phi) */
  );

