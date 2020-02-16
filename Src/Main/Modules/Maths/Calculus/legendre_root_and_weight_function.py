from sympy.solvers import solve
from sympy import legendre, re
from sympy import *

Nmax      = 25 # maximum N
precision = 18 # numerical precision 

x = symbols('x')

print (' /* direct output of this file to legendre_root_and_weight_function.c */\n')
print ('#include "core_lib.h"')
print ('#include "maths_analytic_lib.h"')
print ('#include "maths_general_lib.h"\n')

print ('#define ARRAY_SIZE_LEGENDRE {}'.format(Nmax))
print ('typedef double fdPn_dx_T(const double x);')
print ('fdPn_dx_T *dP_dx[ARRAY_SIZE_LEGENDRE];')
print ('double Legendre_root[ARRAY_SIZE_LEGENDRE][ARRAY_SIZE_LEGENDRE] = {{0}};\n')
# declaration
print ('double Legendre_weight_function(const double x, const unsigned N);')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N);')
print ('double dLegendre_dx(const unsigned n, const double x);')
print ('void init_Legendre_root_function(void);')
print ('void init_dLegendre_dx(void);')
for n in range(Nmax):
  print('static double dP{}_dx(const double x);'.format(n))

# weight function:
print ('\n/* weight = 2./((1-x^2)*(dPn(x)/dx)^2) */')
print ('double Legendre_weight_function(const double x, const unsigned n)')
print ('{')
print ('  if (EQL(x,1.) || EQL(x,-1.))')
print ('    abortEr("Bad argument for Legendre weight function.\\n");')
print ('  return 2./((1-Pow2(x))*Pow2(dLegendre_dx(n,x)));')
print ('}')

# dPn/dx
print ('\n/* ->return value: d(Pn(x))/dx, Pn is legendre(n,x) polynomial */')
print ('double dLegendre_dx(const unsigned n, const double x)')
print ('{')
print ('  if (x > 1. || x < -1.)')
print ('    abortEr("x exceeds from [-1,1] interval.\\n");')
print ('  if (n >= {})'.format(Nmax))
print ('    abortEr("n exceeds the maximum.\\n"\n'\
       '            "To go higher number change Nmax in \'legendre_root_and_weight_function.py\'.\\n");')
print ('  return dP_dx[n](x);')
print ('}')

# roots
print ('\n/* rootN is the n-th root of (legendere(N,x)). the following is like a table.')
print ('// note: rootN starts from 0. */')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N)')
print ('{')
print ('  if (N >= {})'.format(Nmax))
print ('    abortEr("N exceeds from the maximum.\\n");\n')
print ('  if (rootN >= N)') # derivative makes it N-1 and since starts form 0, the max rootN is N-2
print ('    abortEr("root number exceeds from the maximum.\\n");\n')
print ('  return Legendre_root[N][rootN];')
print ('}')


# generating dPn/dx functions
for n in range(Nmax):
    print('\n/* dP{}/dx */'.format(n))
    print('static double dP{}_dx(const double x)'.format(n))
    print('{')
    if n == 0 or n == 1:
        print ('  UNUSED(x);')
    print('  return '+ ccode(simplify(diff(legendre(n, x),x))) + ';')
    print('}\n')


# dPn/dx
print ('\n/* initializing d(Pn(x))/dx, Pn is legendre(n,x) polynomial */')
print ('void init_dLegendre_dx(void)')
print ('{')
# assign functions
for n in range(Nmax):
  print ('  dP_dx[{0}] = dP{0}_dx;'.format(n))        
print ('}')

# roots
print ('/* initializing root table */')
print ('void init_Legendre_root_function(void)')
print ('{')
print ('  Legendre_root[0][0] = 0.;')
for i in range(1,Nmax):
    root  = solve(legendre(i,x),x,minimal=True,quick=True,warn=True)
    Realroot  = []
    Nroot = len(root)
    for j in range(Nroot):# make it real, some of the roots has 10-17 I part
      Realroot.append(re(N(root[j],precision)));
    Realroot  = sorted(Realroot)
    for j in range(Nroot):
      print ('  Legendre_root[{}][{}] = {};'.format(i,j,Realroot[j]))
print ('}')

