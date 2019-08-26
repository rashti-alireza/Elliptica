from sympy import assoc_legendre
from sympy import *

lmax = 15 # maximum l for P_{l}^{m}

x=symbols("x")

print (' /* direct output of this file to associated_legendre.c */')

print ('#include "core_lib.h"')
print ('#include "maths_analytic_lib.h"\n')
print ('#define LMAX_ARRAY_SIZE_PLM {}'.format(lmax))
print ('typedef double fPlm_T (const double x);')
print ('fPlm_T *P[LMAX_ARRAY_SIZE_PLM][LMAX_ARRAY_SIZE_PLM];')
print ('fPlm_T *P_[LMAX_ARRAY_SIZE_PLM][LMAX_ARRAY_SIZE_PLM];') # for negative m, this is more accurate 
                                    # as oppose to using relation between P_{l}^{m} and P_{l}^{-m}
# declaration
print ('void init_associated_legendre(void);')
print ('double associated_legendre(const int l, const int m, const double x);')
for l in range(lmax):
    for m in range(l+1):
        print('static double associated_legendre_P_l{}m{}(const double x);'.format(l,m))

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('static double associated_legendre_P_l{}m_{}(const double x);'.format(l,-m))
        
print ('\n\n')        
print ('/* initializing P_{l}^{m}(x) table */')
print ('void init_associated_legendre(void)')
print ('{')
# assign functions
for l in range(lmax):
    for m in range(l+1):
        print ('  P[{0}][{1}] = associated_legendre_P_l{0}m{1};'.format(l,m))
        
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print ('  P_[{0}][{1}] = associated_legendre_P_l{0}m_{1};'.format(l,-m))        
print ('}')

print ('/* P_{l}^{m}(x)=\\left( -1\\right) ^{m}\\left( 1-x^{2}\\right) ^{\\frac {m} {2}}\\frac {d^{m}P_{l}\\left( x\\right) } {dx^{m}}')
print ('// ->return value: P_{l}^{m}(x) */')
print ('double associated_legendre(const int l, const int m, const double x)')
print ('{')
print ('  if (x > 1. || x < -1.)')
print ('    abortEr("x exceeds from [-1,1] interval.\\n");')
print ('  if (l >= {})'.format(lmax))
print ('    abortEr("l exceeds the maximum.\\n"\n'\
       '            "To go higher number change lmax in \'associated_legendre.py\'.\\n");')
print ('  if (l < 0)')
print ('    abortEr("l is negative.\\n");')
print ('  if (abs(m) > l)')
print ('    return 0;')
print ('  if (m < 0)')
print ('    return P_[l][-m](x);')
print ('  return P[l][m](x);')
print ('}')


for l in range(lmax):
    for m in range(l+1):
        print('\n/* P_{{{}}}^{{{}}} */'.format(l,m))
        print('static double associated_legendre_P_l{}m{}(const double x)'.format(l,m))
        print('{')
        if l == 0 and m == 0:
            print ('  UNUSED(x);')
        print('  return '+ ccode(assoc_legendre(l,m, x)) + ';')
        print('}\n')

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('\n/* P_{{{}}}^{{{}}} */'.format(l,m))
        print('static double associated_legendre_P_l{}m_{}(const double x)'.format(l,-m))
        print('{')
        print('  return '+ ccode(assoc_legendre(l,m, x)) + ';')
        print('}\n')
        