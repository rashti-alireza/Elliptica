### Calculating free data _gamma,_gammaI,_dgamma,_ddgamma,_dddgamma
## analytically.

#########################
### H O W  T O  U S E ###
#########################
# Issue the following (TAKES A VERY VERY LONG TIME):
# 
# $ cpi bbn_free_data_g_analytic.math > bbn_free_data_g_analytic.c && \
#   sed -i '/Welcome to Cpi/,$d' bbn_free_data_g_analytic.c


`from __future__ import division`
`from sympy import *`
`from sympy.tensor.tensor import TensorIndexType, TensorHead`
`from sympy.tensor.tensor import tensor_indices`
`from sympy.functions import transpose`
`import sys`
`import re`
`import os`

### ccode generator
`def mcode(m):`
`  code = ""`
`  if (0):`
`    pid = os.getpid()`
`    mfile_name=".mfile_temp_{}".format(pid)`
`    cfile_name=".cfile_temp_{}".format(pid)`
`    rhs  =""`
`    try:`
`      rhs  = mathematica_code((m).doit())`
`    except:`
`      rhs  = mathematica_code(m)`
``
`    math_code  = rhs + ";\n"`
`    math_code += "Simplify[%,TimeConstraint->2000];\n"`
`    math_code += "c=CForm[%];\n"`
`    math_code += "Print[c]"+"\n"`
## write into mathematica file
`    mfile = open(mfile_name,"w")`
`    mfile.write(math_code)`
`    mfile.close()`

## run mathematica and output into another
`    cmd = "math -run -noprompt < {} 1> {}".format(mfile_name,cfile_name)`
`    ret=os.system(cmd)`

## read math results.
`    cfile=open(cfile_name,"r")`
`    code =cfile.read()`
`    cfile.close()`

## delete the files
`    cmd = "rm -rf {} {}".format(mfile_name,cfile_name)`
`    os.system(cmd)`
`  else:`
`    code = ccode(m)`
``
## return the code and break long lines
`  code = '{};'.format(code)`
`  code = re.sub(r'\b1.0\*\b','',code)`
`  code = re.sub(r'(.{50,70}[\+\-\*/]+\s?)','\\1{}n'.format("\\"),code)`
`  sys.stdout.flush()`
`  return code`
``
``
`### coords and parameters`
`x,y,z = symbols('x,y,z')`
``
`r0,M_BH,a_BH,Lambda = symbols('r0,M_BH,a_BH,Lambda')`
``
`BH_center_x,BH_center_y,BH_center_z = \`
`  symbols('BH_center_x,BH_center_y,BH_center_z')`
`  `
`chi,chi_U0,chi_U1,chi_U2 = symbols('chi,chi_U0,chi_U1,chi_U2')`
`x_CM,y_CM,Omega_BHNS = symbols('x_CM,y_CM,Omega_BHNS')`
``
`Bx,By,Bz,B2 = symbols('Bx,By,Bz,B2')`
`phiy,phiz   = symbols('phiy,phiz')`
``
`r = (x**2+y**2+z**2)**0.5`
`e = exp(-(r/r0)**4)`
``
`### functions`
`rKS = symbols('rKS',cls=Function)`
`HKS = symbols('HKS',cls=Function)`
`TBR = symbols('TBR',cls=Function)`
`PB  = symbols('TB',cls=Function)`
`PR  = symbols('TR',cls=Function)`
``
``
`### r funciont in Kerr-Schild coords:`
`def rKS(_x,_y,_z,a):`
`  r2 = _x**2+_y**2+_z**2`
`  a2 = a**2`
`  return (0.5*(r2-a2+sqrt((r2-a2)**2+4*a2*_z**2)))`
``
``
`### H function in Kerr-Schild coords`
`def HKS (M_BH,_r,a,_z,Lambda):`
`  _k2 = _z/_r;`
`  a2  = a**2`
`  _r2 = _r**2`
`  return Lambda*M_BH*_r/(_r2+a2*_k2**2)`
``
``
`### boost and rotation transformation, dir = [1,-1]`
`def TBR(B,Ry,Rz,dir,v):`
`  if (dir == 1):`
`    return (B*Rz*Ry*v)`
`  elif (dir == -1):`
`    return (Ry**-1 * Rz**-1 *B**-1 *v)`
`  else:`
`    raise Exception('Bad argument')`
``
`  `
``
`### populate boost transformation:`
`def PB(Bx,By,Bz,B2):`
`   gamma  = (1-B2)**(-0.5)`
`   b = Matrix([`
`     [gamma    , -gamma*Bx              ,-gamma*By              ,-gamma*Bz],`
`     [-gamma*Bx, 1+(gamma-1)*Bx**2/B2,(gamma-1)*Bx*By/B2     ,(gamma-1)*Bx*Bz/B2],`
`     [-gamma*By, (gamma-1)*Bx*By/B2     ,1+(gamma-1)*By**2/B2,(gamma-1)*By*Bz/B2],`
`     [-gamma*Bz, (gamma-1)*Bx*Bz/B2     ,(gamma-1)*By*Bz/B2     ,1+(gamma-1)*Bz**2/B2]`
`     ])`
`     `
`   return b;`
``
`### populate rotation transformaion`
`def PR(phiy,phiz):`
`  cy = cos(phiy)`
`  sy = sin(phiy)`
`  cz = cos(phiz)`
`  sz = sin(phiz)`
`  `
`  ty = Matrix([`
`     [1,0 ,0,0],`
`     [0,cy,0,sy],`
`     [0,0 ,1,0],`
`     [0,-sy,0,cy]`
`     ])`
`     `
`  tz = Matrix([`
`     [1,0 ,0  ,0],`
`     [0,cz,-sz,0],`
`     [0,sz,cz ,0],`
`     [0,0 ,0  ,1]`
`     ])`
`     `
`  return ty, tz  `
``
``
`### get ks and H needed to populate Kerr-Schild`
`def get_ks_and_H(x,y,z,a,m,Bx,By,Bz,B2,phiy,phiz):`
`  x_mu     = transpose(Matrix([[0,x,y,z]]))`
`  tB       = PB(Bx,By,Bz,B2)`
`  tRy, tRz = PR(phiy,phiz)`
`  _x_mu = TBR(tB,tRy,tRz,-1,x_mu)`
`  _x    = simplify(_x_mu[1,0])`
`  _y    = simplify(_x_mu[2,0])`
`  _z    = simplify(_x_mu[3,0])`
`  _r    = rKS(_x,_y,_z,a)`
`  a2  = a**2`
`  _r2 = _r**2`
`  _k0 = simplify((_r*_x+a*_y)/(_r2+a2))`
`  _k1 = simplify((_r*_y-a*_x)/(_r2+a2))`
`  _k2 = simplify(_z/_r)`
`  _kt = 1`
`  `
`  _k_mu = transpose(Matrix([[_kt,_k0,_k1,_k2]]))`
`  k_mu  = TBR(tB,tRy,tRz,1,_k_mu)`
`  `
`  print("/* kt */")`
`  sys.stdout.flush()`
`  kt = simplify(k_mu[0,0])`
`  `
`  print("/* k0 */")`
`  sys.stdout.flush()`
`  k0 = simplify(k_mu[1,0])`
`  `
`  print("/* k1 */")`
`  sys.stdout.flush()`
`  k1 = simplify(k_mu[2,0])`
`  `
`  print("/* k2 */")`
`  sys.stdout.flush()`
`  k2 = simplify(k_mu[3,0])`
`  `
`  print("/* H */")`
`  sys.stdout.flush()`
`  H = simplify(HKS(m,_r,a,_z,Lambda))`
`  `
`  return kt,k0,k1,k2,H`
``
`### evaluates:`
`kt,k0,k1,k2,H = get_ks_and_H(x,y,z,a_BH,M_BH,Bx,By,Bz,B2,phiy,phiz)`
`r = (x**2+y**2+z**2)**0.5`
`e = exp(-(r/r0)**4)`
`C = 2*H*e`
`print("/* A */")`
`sys.stdout.flush()`
`A = simplify(1/(1+C*(k0**2+k1**2+K2**2)))`
``
`_g00 = 1.+C*k0*k0`
`_g01 = C*k0*k1`
`_g02 = C*k0*k2`
`_g11 = 1+C*k1*k1`
`_g12 = C*k1*k2`
`_g22 = 1+C*k2*k2`
``
`_g = Matrix([`
`      [_g00,_g01,_g02],`
`      [_g01,_g11,_g12],`
`      [_g02,_g12,_g22]`
`      ])`
`_ig = inverse(_g)`
``
`########################`
`### C code generator ###`
`########################`
``
print('#include "bbn_headers.h"')
print('void bbn_free_data_g_gI_analytic(Patch_T *const patch);')
print('void bbn_free_data_g_gI_analytic(Patch_T *const patch)')
print('{')
print('    unsigned nn,ijk;')
print('    nn = patch->nn;')
print('\n')
print('    REALLOC_v_WRITE_v(_gamma_D2D2)')
print('    REALLOC_v_WRITE_v(_gamma_D0D2)')
print('    REALLOC_v_WRITE_v(_gamma_D0D0)')
print('    REALLOC_v_WRITE_v(_gamma_D0D1)')
print('    REALLOC_v_WRITE_v(_gamma_D1D2)')
print('    REALLOC_v_WRITE_v(_gamma_D1D1)')
print('    REALLOC_v_WRITE_v(_gammaI_U0U2)')
print('    REALLOC_v_WRITE_v(_gammaI_U0U0)')
print('    REALLOC_v_WRITE_v(_gammaI_U0U1)')
print('    REALLOC_v_WRITE_v(_gammaI_U1U2)')
print('    REALLOC_v_WRITE_v(_gammaI_U1U1)')
print('    REALLOC_v_WRITE_v(_gammaI_U2U2)')
print('    ')
print('    for (ijk = 0; ijk < nn; ++ijk)')
print('    {')
print('      double x,y,z,r,H,k0,k1,k2,kt;')
print('\n')
print('      x = patch->node[ijk]->x[0]-BH_center_x;')
print('      y = patch->node[ijk]->x[1]-BH_center_y;')
print('      z = patch->node[ijk]->x[2]-BH_center_z;')
print('      r = sqrt(Pow2(x)+Pow2(y)+Pow2(z));')
print('\n')
print('      _gamma_D0D0[ijk] =')
print(mcode(_g[0,0])

print('      _gamma_D0D1[ijk] =') 
print(mcode(_g[0,1])

print('      _gamma_D0D2[ijk] =')
print(mcode(_g[0,2])

print('      _gamma_D1D1[ijk] =')
print(mcode(_g[1,1])

print('      _gamma_D1D2[ijk] =')
print(mcode(_g[1,2])

print('      _gamma_D2D2[ijk] =')
print(mcode(_g[2,2])

print('\n')
print('      _gammaI_U0U0[ijk] =')
print(mcode(_ig[0,0])

print('      _gammaI_U0U1[ijk] =')
print(mcode(_ig[0,1])

print('      _gammaI_U0U2[ijk] =')
print(mcode(_ig[0,2])

print('      _gammaI_U1U1[ijk] =')
print(mcode(_ig[1,1])

print('      _gammaI_U1U2[ijk] =')
print(mcode(_ig[1,2])

print('      _gammaI_U2U2[ijk] =')
print(mcode(_ig[2,2])

print('      /* quick test check _gamma * _gammaI = delta */')
print('      if (0)')
print('      {')
print('          double delta_U0D0 = ')
print('        _gammaI_U0U0[ijk]*_gamma_D0D0[ijk] + _gammaI_U0U1[ijk]*')
print('        _gamma_D0D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D0D2[ijk];')
print('\n')
print('          double delta_U0D1 = ')
print('        _gammaI_U0U0[ijk]*_gamma_D0D1[ijk] + _gammaI_U0U1[ijk]*')
print('        _gamma_D1D1[ijk] + _gammaI_U0U2[ijk]*_gamma_D1D2[ijk];')
print('\n')
print('          double delta_U0D2 = ')
print('        _gammaI_U0U0[ijk]*_gamma_D0D2[ijk] + _gammaI_U0U1[ijk]*')
print('        _gamma_D1D2[ijk] + _gammaI_U0U2[ijk]*_gamma_D2D2[ijk];')
print('\n')
print('          double delta_U1D2 = ')
print('        _gammaI_U0U1[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U1[ijk]*')
print('        _gamma_D1D2[ijk] + _gammaI_U1U2[ijk]*_gamma_D2D2[ijk];')
print('\n')
print('          double delta_U1D0 = ')
print('        _gammaI_U0U1[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U1[ijk]*')
print('        _gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D0D2[ijk];')
print('\n')
print('         double delta_U1D1 = ')
print('        _gammaI_U0U1[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U1[ijk]*')
print('        _gamma_D1D1[ijk] + _gammaI_U1U2[ijk]*_gamma_D1D2[ijk];')
print('        ')
print('          double delta_U2D2 = ')
print('        _gammaI_U0U2[ijk]*_gamma_D0D2[ijk] + _gammaI_U1U2[ijk]*')
print('        _gamma_D1D2[ijk] + _gammaI_U2U2[ijk]*_gamma_D2D2[ijk];')
print('\n')
print('          double delta_U2D0 = ')
print('        _gammaI_U0U2[ijk]*_gamma_D0D0[ijk] + _gammaI_U1U2[ijk]*')
print('        _gamma_D0D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D0D2[ijk];')
print('\n')
print('          double delta_U2D1 = ')
print('        _gammaI_U0U2[ijk]*_gamma_D0D1[ijk] + _gammaI_U1U2[ijk]*')
print('        _gamma_D1D1[ijk] + _gammaI_U2U2[ijk]*_gamma_D1D2[ijk];')
print('\n')
print('        if(!EQL(delta_U1D1,1)||!isfinite(delta_U1D1))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U0D1,0)||!isfinite(delta_U0D1))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U0D2,0)||!isfinite(delta_U0D2))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U1D2,0)||!isfinite(delta_U1D2))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U0D0,1)||!isfinite(delta_U0D0))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U2D1,0)||!isfinite(delta_U2D1))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U2D2,1)||!isfinite(delta_U2D2))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U2D0,0)||!isfinite(delta_U2D0))  Error0("_gammaI is not correct!\n");')
print('        if(!EQL(delta_U1D0,0)||!isfinite(delta_U1D0))  Error0("_gammaI is not correct!\n");')
print('\n')
print('      }')
print('    }')
print('}')
