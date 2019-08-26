from sympy import Ynm
from sympy import *
import re

lmax = 15 # maximum l for Y_{l}^{m}

phi=symbols("phi")
theta=symbols("theta")

print (' /* direct output of this file to spherical_harmonic_Ylm.c */')

print ('#include "core_lib.h"')
print ('#include "maths_analytic_lib.h"')
print ('#include "maths_general_lib.h"')
print ('#include <complex.h>\n')
print ('static const int lmax = {};'.format(lmax))
print ('typedef double complex fYlm_T (const double theta, const double phi);\n')
# declaration
print ('double complex Ylm(const unsigned l, const int m, const double theta, const double phi);')
for l in range(lmax):
    for m in range(l+1):
        print('double complex Ylm_l{}m{}(const double theta, const double phi);'.format(l,m))

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('double complex Ylm_l{}m_{}(const double theta, const double phi);'.format(l,-m))
        
print ('\n\n')        
print ('/* Y_n^m(\\theta, \\varphi) := \\sqrt{\\frac{(2n+1)(n-m)!}{4\\pi(n+m)!}} exp(i m \\varphi)\\mathrm{P}_n^m\\left(\\cos(\\theta)\\right) ');
print ('// ->return value: Y_{l}^{m}(x) */')
print ('double complex Ylm(const unsigned l, int m, const double theta, const double phi)')
print ('{')
print ('  fYlm_T *Y[lmax][lmax];')
print ('  fYlm_T *Y_[lmax][lmax];') # for negative m, this is more accurate 
                                    # as oppose to using relation between P_{l}^{m} and P_{l}^{-m}
print ('  if (theta > M_PI || theta < 0)')
print ('    abortEr("theta exceeds from [0,pi] interval.\\n");')
print ('  if (phi > 2*M_PI || phi < 0)')
print ('    abortEr("phi exceeds from [0,2*pi] interval.\\n");')
print ('  if (l >= {})'.format(lmax))
print ('    abortEr("l exceeds the maximum.\\n");')
print ('  if ((unsigned)abs(m) > l)')
print ('    return 0;')

# assign functions
for l in range(lmax):
    for m in range(l+1):
        print ('  Y[{0}][{1}] = Ylm_l{0}m{1};'.format(l,m))
        
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print ('  Y_[{0}][{1}] = Ylm_l{0}m_{1};'.format(l,-m))        

print ('  if (m < 0)')
print ('    return Y_[l][-m](theta,phi);')
print ('  return Y[l][m](theta,phi);')
print ('}')

for l in range(lmax):
    for m in range(0,l+1,1):
        print('\n/* Y_{{{}}}^{{{}}} */'.format(l,m))
        print('double complex Ylm_l{}m{}(const double theta, const double phi)'.format(l,m))
        print('{')
        if l == 0:
            print ('  UNUSED(theta);')
        if m == 0:
            print ('  UNUSED(phi);')
        expr = ccode(simplify(Ynm(l, m, theta, phi).expand(func=True)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('\n/* Y_{{{}}}^{{{}}} */'.format(l,m))
        print('double complex Ylm_l{}m_{}(const double theta, const double phi)'.format(l,-m))
        print('{')
        expr = ccode(simplify(Ynm(l, m, theta, phi).expand(func=True)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')
        