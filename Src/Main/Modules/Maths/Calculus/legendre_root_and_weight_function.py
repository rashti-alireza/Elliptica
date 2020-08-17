from sympy.solvers import solve
from sympy import legendre, re
from sympy import *

Nmax      = 100 # maximum N
precision = 18 # numerical precision 

x = symbols('x')

print (' /* direct output of this file to legendre_root_and_weight_function.c */\n')
print ('#include "core_lib.h"')
print ('#include "maths_special_functions_lib.h"')
print ('#include "maths_general_lib.h"\n')
print ('#define ARRAY_SIZE_LEGENDRE {}'.format(Nmax))
print ('double Legendre_root[ARRAY_SIZE_LEGENDRE][ARRAY_SIZE_LEGENDRE] = {{0}};\n')
# declaration
print ('double Legendre_weight_function(const double x, const unsigned N);')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N);')
print ('void init_Legendre_root_function(void);')

# weight function:
print ('\n/* weight = 2./(dPn(cos(theta))/dtheta)^2) */')
print ('double Legendre_weight_function(const double x, const unsigned n)')
print ('{')
print ('  const double dp_dth = d_associated_legendre_dtheta((int)n,0,x);')
print ('  return 2./(Pow2(dp_dth));')
print ('}')

# roots
print ('\n/* rootN is the n-th root of (legendere(N,x)). the following is like a table.')
print ('// note: rootN starts from 0. */')
print ('double Legendre_root_function(const unsigned rootN, const unsigned N)')
print ('{')
print ('  if (N >= {})'.format(Nmax))
print ('    Error0("N exceeds from the maximum.\\n");\n')
print ('  if (rootN >= N)') # derivative makes it N-1 and since starts form 0, the max rootN is N-2
print ('    Error0("root number exceeds from the maximum.\\n");\n')
print ('  return Legendre_root[N][rootN];')
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

