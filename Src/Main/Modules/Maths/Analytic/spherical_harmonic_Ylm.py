from sympy import Ynm
from sympy import *
import re

lmax = 15 # maximum l for Y_{l}^{m}

phi   = symbols("phi")
theta = symbols("theta")

print (' /* direct output of this file to spherical_harmonic_Ylm.c */')

print ('#include "core_lib.h"')
print ('#include "maths_analytic_lib.h"')
print ('#include "maths_general_lib.h"')
print ('#include <complex.h>\n')
print ('#define LMAX_ARRAY_SIZE_YLM {}'.format(lmax))
print ('typedef double complex fYlm_T (const double theta, const double phi);')
print ('fYlm_T *Y[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];')
print ('fYlm_T *Y_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];') # for negative m, this is more accurate 
                                  # as oppose to using relation between P_{l}^{m} and P_{l}^{-m}
                                  
print ('fYlm_T *dY_dphi[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];')
print ('fYlm_T *dY_dphi_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];') # for negative m, this is more accurate 
                                        # as oppose to using relation between P_{l}^{m} and P_{l}^{-m}
                                          
print ('fYlm_T *dY_dtheta[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];')
print ('fYlm_T *dY_dtheta_[LMAX_ARRAY_SIZE_YLM][LMAX_ARRAY_SIZE_YLM];\n') # for negative m, this is more accurate 
                                        # as oppose to using relation between P_{l}^{m} and P_{l}^{-m}
# declarations:
print ('void init_Ylm(void);')
print ('void init_dYlm_dphi(void);')
print ('void init_dYlm_dtheta(void);')
print ('double complex Ylm(const unsigned l, const int m, const double theta, const double phi);')
print ('double complex dYlm_dphi(const unsigned l, const int m, const double theta, const double phi);')
print ('double complex dYlm_dtheta(const unsigned l, const int m, const double theta, const double phi);')
# Ylm:
for l in range(lmax):
    for m in range(l+1):
        print('static double complex Ylm_l{}m{}(const double theta, const double phi);'.format(l,m))

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('static double complex Ylm_l{}m_{}(const double theta, const double phi);'.format(l,-m))
# dYlm:
for l in range(lmax):
    for m in range(l+1):
        print('static double complex dYlm_dphi_l{}m{}(const double theta, const double phi);'.format(l,m))

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('static double complex dYlm_dphi_l{}m_{}(const double theta, const double phi);'.format(l,-m))

# dYlm:
for l in range(lmax):
    for m in range(l+1):
        print('static double complex dYlm_dtheta_l{}m{}(const double theta, const double phi);'.format(l,m))

for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('static double complex dYlm_dtheta_l{}m_{}(const double theta, const double phi);'.format(l,-m))

# Ylm:        
print ('\n\n')
print ('/* initializing Y_{l}^{m}(x) */')
print ('void init_Ylm(void)')
print ('{')

# assign functions
for l in range(lmax):
    for m in range(l+1):
        print ('  Y[{0}][{1}] = Ylm_l{0}m{1};'.format(l,m))
        
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print ('  Y_[{0}][{1}] = Ylm_l{0}m_{1};'.format(l,-m))        

print ('}')


print ('/* Y_n^m(\\theta, \\varphi) := \\sqrt{\\frac{(2n+1)(n-m)!}{4\\pi(n+m)!}} exp(i m \\varphi)\\mathrm{P}_n^m\\left(\\cos(\\theta)\\right) ');
print ('// ->return value: Y_{l}^{m}(x) */')
print ('double complex Ylm(const unsigned l, int m, const double theta, const double phi)')
print ('{')
print ('  if (theta > M_PI || theta < 0)')
print ('    abortEr("theta exceeds from [0,pi] interval.\\n");')
print ('  if (phi > 2*M_PI || phi < 0)')
print ('    abortEr("phi exceeds from [0,2*pi] interval.\\n");')
print ('  if (l >= {})'.format(lmax))
print ('    abortEr("l exceeds the maximum.\\n");')
print ('  if ((unsigned)abs(m) > l)')
print ('    return 0;')
print ('  if (m < 0)')
print ('    return Y_[l][-m](theta,phi);')
print ('  return Y[l][m](theta,phi);')
print ('}')

# dYlm/dphi:
print ('\n\n')        
print ('/* initializing dYlm_dphi */')
print ('void init_dYlm_dphi(void)')
print ('{')
# assign functions
for l in range(lmax):
    for m in range(l+1):
        print ('  dY_dphi[{0}][{1}] = dYlm_dphi_l{0}m{1};'.format(l,m))
        
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print ('  dY_dphi_[{0}][{1}] = dYlm_dphi_l{0}m_{1};'.format(l,-m))        

print ('}')

# dYlm/dphi:
print ('/* ->return value: dY_{l}^{m}(x)/dphi */')
print ('double complex dYlm_dphi(const unsigned l, int m, const double theta, const double phi)')
print ('{')
print ('  if (theta > M_PI || theta < 0)')
print ('    abortEr("theta exceeds from [0,pi] interval.\\n");')
print ('  if (phi > 2*M_PI || phi < 0)')
print ('    abortEr("phi exceeds from [0,2*pi] interval.\\n");')
print ('  if (l >= {})'.format(lmax))
print ('    abortEr("l exceeds the maximum.\\n");')
print ('  if ((unsigned)abs(m) > l)')
print ('    return 0;')
print ('  if (m < 0)')
print ('    return dY_dphi_[l][-m](theta,phi);')
print ('  return dY_dphi[l][m](theta,phi);')
print ('}')

# dYlm/dtheta:
print ('\n\n')        
print ('/* initialize dYlm_dtheta */')
print ('void init_dYlm_dtheta(void)')
print ('{')
# assign functions
for l in range(lmax):
    for m in range(l+1):
        print ('  dY_dtheta[{0}][{1}] = dYlm_dtheta_l{0}m{1};'.format(l,m))
        
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print ('  dY_dtheta_[{0}][{1}] = dYlm_dtheta_l{0}m_{1};'.format(l,-m))        

print ('}')
# dYlm/theta:
print ('/* ->return value: dY_{l}^{m}(x)/dtheta */')
print ('double complex dYlm_dtheta(const unsigned l, int m, const double theta, const double phi)')
print ('{')
print ('  if (theta > M_PI || theta < 0)')
print ('    abortEr("theta exceeds from [0,pi] interval.\\n");')
print ('  if (phi > 2*M_PI || phi < 0)')
print ('    abortEr("phi exceeds from [0,2*pi] interval.\\n");')
print ('  if (l >= {})'.format(lmax))
print ('    abortEr("l exceeds the maximum.\\n");')
print ('  if ((unsigned)abs(m) > l)')
print ('    return 0;')
print ('  if (m < 0)')
print ('    return dY_dtheta_[l][-m](theta,phi);')
print ('  return dY_dtheta[l][m](theta,phi);')
print ('}')

# Y_{l}^{m}, m >= 0
for l in range(lmax):
    for m in range(0,l+1,1):
        print('\n/* Y_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex Ylm_l{}m{}(const double theta, const double phi)'.format(l,m))
        print('{')
        if l == 0:
            print ('  UNUSED(theta);')
        if m == 0:
            print ('  UNUSED(phi);')
        expr = ccode(simplify(Ynm(l, m, theta, phi).expand(func=True)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

# Y_{l}^{m}, m < 0
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('\n/* Y_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex Ylm_l{}m_{}(const double theta, const double phi)'.format(l,-m))
        print('{')
        expr = ccode(simplify(Ynm(l, m, theta, phi).expand(func=True)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

# dY_dphi_{l}^{m}, m >= 0
for l in range(lmax):
    for m in range(0,l+1,1):
        print('\n/* dY_dphi_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex dYlm_dphi_l{}m{}(const double theta, const double phi)'.format(l,m))
        print('{')
        if m == 0: # derivative becomes 0
            print ('  UNUSED(phi);')
            print ('  UNUSED(theta);')
        expr = ccode(simplify(diff(Ynm(l, m, theta, phi).expand(func=True),phi)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

# dY_dphi_{l}^{m}, m < 0
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('\n/* dY_dphi_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex dYlm_dphi_l{}m_{}(const double theta, const double phi)'.format(l,-m))
        print('{')
        expr = ccode(simplify(diff(Ynm(l, m, theta, phi).expand(func=True),phi)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

# dY_dtheta_{l}^{m}, m >= 0
for l in range(lmax):
    for m in range(0,l+1,1):
        print('\n/* dY_dtheta_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex dYlm_dtheta_l{}m{}(const double theta, const double phi)'.format(l,m))
        print('{')
        if l == 0:
            print ('  UNUSED(phi);')
            print ('  UNUSED(theta);')
        elif m == 0:
            print ('  UNUSED(phi);')
        expr = ccode(simplify(diff(Ynm(l, m, theta, phi).expand(func=True),theta)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')

# dY_dtheta_{l}^{m}, m < 0
for l in range(lmax):
    for m in range(-1,-l-1,-1):
        print('\n/* dY_dtheta_{{{}}}^{{{}}} */'.format(l,m))
        print('static double complex dYlm_dtheta_l{}m_{}(const double theta, const double phi)'.format(l,-m))
        print('{')
        expr = ccode(simplify(diff(Ynm(l, m, theta, phi).expand(func=True),theta)))
        expr = re.sub(r'\bexp\(','cexp(',expr);
        print('  return '+ expr + ';')
        print('}\n')
        