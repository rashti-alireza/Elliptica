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

print ('static const int Nmax = {};'.format(Nmax))
print ('typedef double fdPn_dx_T(const double x);\n')

# declaration
print ('double Legendre_weight_function(const double x, const unsigned N);')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N);')
print ('double dlegendre_dx(const unsigned n, const double x);')
for n in range(Nmax):
  print('static double dP{}_dx(const double x);'.format(n))
  

# weight function:
print ('\n/* weight = 2./((1-x^2)*(dPn(x)/dx)^2) */')
print ('double Legendre_weight_function(const double x, const unsigned n)')
print ('{')
print ('  if (EQL(x,1.) || EQL(x,-1.))')
print ('    abortEr("Bad argument for Legendre weight function.\\n");')
print ('  return 2./((1-SQR(x))*SQR(dlegendre_dx(n,x)));')
print ('}')

# dPn/dx
print ('\n/* ->return value: d(Pn(x))/dx, Pn is legendre(n,x) polynomial */')
print ('double dlegendre_dx(const unsigned n, const double x)')
print ('{')
print ('  fdPn_dx_T *dP_dx[Nmax];')
print ('  if (x > 1. || x < -1.)')
print ('    abortEr("x exceeds from [-1,1] interval.\\n");')
print ('  if (n >= {})'.format(Nmax))
print ('    abortEr("n exceeds the maximum.\\n"\n'\
       '            "To go higher number change Nmax in \'legendre_root_and_weight_function.py\'.\\n");')
# assign functions
for n in range(Nmax):
  print ('  dP_dx[{0}] = dP{0}_dx;'.format(n))        
print ('  return dP_dx[n](x);')
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

# roots
print ('\n/* rootN is the n-th root of (legendere(N,x)). the following is like a table.')
print ('// note: rootN starts from 0. */')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N)')
print ('{')
print ('  double root[{}][{}] = {{0}};'.format(Nmax,Nmax))
print ('  if (N >= {})'.format(Nmax))
print ('    abortEr("N exceeds from the maximum.\\n");\n')
print ('  if (rootN >= N)') # derivative makes it N-1 and since starts form 0, the max rootN is N-2
print ('    abortEr("root number exceeds from the maximum.\\n");\n')
print ('  root[0][0] = 0.;')

for i in range(1,Nmax):
    root  = solve(legendre(i,x),x,minimal=True,quick=True,warn=True)
    Nroot = len(root)
    for j in range(Nroot):
      print ('  root[{}][{}] = {};'.format(i,j,re(N(root[j],precision))))

print ('  return root[N][rootN];')
print ('}')

