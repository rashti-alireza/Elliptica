from sympy.solvers import solve
from sympy import legendre, re
from sympy import *

Nmax      = 35 # maximum N
precision = 18 # numerical precision 

x = symbols('x')


print (' /* direct output of this file to lobatto_root_and_weight_function.c */\n')
print ('#include "core_lib.h"')
print ('#include "maths_analytic_lib.h"')
print ('#include "maths_general_lib.h"\n')
print ('double Lobbatto_weight_function(const double x, const unsigned N);')
print ('double Lobbatto_root_function(const unsigned rootN, const unsigned N);\n')

print ('/* weight = 2/(n*(n-1)*(P_{n-1}(x))^2), if x == +/- 1 w = 2/(n*(n-1)) */')
print ('double Lobbatto_weight_function(const double x, const unsigned n)')
print ('{')
print ('  double fac = 2./(n*(n-1));')
print ('  if (EQL(x,1.) || EQL(x,-1.))')
print ('    return fac;')
print ('  return 2./(n*(n-1)*SQR(associated_legendre((int)n-1,0,x)));')
print ('}')

print ('/* rootN is the n-th root of diff(legendere(N,x)). the following is like a table.')
print ('// note: rootN starts from 0. */')
print ('double Lobbatto_root_function(const unsigned rootN, const unsigned N)')
print ('{')
print ('  double root[{}][{}] = {{0}};'.format(Nmax,Nmax))
print ('  if (N >= {})'.format(Nmax))
print ('    abortEr("N exceeds from the maximum.\\n");\n')
print ('  if (rootN >= N-1)') # derivative makes it N-1 and since starts form 0, the max rootN is N-2
print ('    abortEr("root number exceeds from the maximum.\\n");\n')
print ('  root[0][0] = 0.;')
print ('  root[1][0] = 0.;')
print ('  root[2][0] = 0.;')

for i in range(3,Nmax):
    root  = solve(diff(legendre(i,x),x),x,minimal=True,quick=True,warn=True)
    Nroot = len(root)
    for j in range(Nroot):
      print ('  root[{}][{}] = {};'.format(i,j,re(N(root[j],precision))))

print ('  return root[N][rootN];')
print ('}')

