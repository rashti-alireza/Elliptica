/*
// Alireza Rashti
// July 2018
*/

#include "analytic_functions.h"

/* the followings are some analytic functions used for various purposes.
// ->return value: a pointer to value of function on whole grid.
*/

/* f: grid -> x */
double *x_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    /* go over all patch's points */
    for (i = in; i < fi; ++i)
      f[i] = x_(i);
      
    in = fi;
  }
  return f;
}

/* f: grid -> y */
double *y_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = y_(i);
      
    in = fi;
  }
  return f;
}

/* f: grid -> z */
double *z_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = z_(i);
      
    in = fi;
  }
  return f;
}

/* f: grid -> r = sqrt(x^2+y^2+z^2) */
double *r_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sqrt(pow(x_(i),2)+pow(y_(i),2)+pow(z_(i),2));
      
    in = fi;
  }
  return f;
}

/* f: grid -> cos(x) */
double *cosx_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = cos(x_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> cos(y) */
double *cosy_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = cos(y_(i));
      
    in = fi;
  }
  return f;
}
/* f: grid -> cos(z) */
double *cosz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = cos(z_(i));
      
    in = fi;
  }
  return f;
}
/* f: grid -> sin(x) */
double *sinx_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sin(x_(i));
      
    in = fi;
  }
  return f;
}
/* f: grid -> sin(y) */
double *siny_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sin(y_(i));
      
    in = fi;
  }
  return f;
}


/* f: grid -> sin(z) */
double *sinz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sin(z_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> cos(xyz) */
double *cosxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = cos(x_(i)*y_(i)*z_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> cos^4(xyz) */
double *cos4xyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = pow(cos(x_(i)*y_(i)*z_(i)),4);
      
    in = fi;
  }
  return f;
}

/* f: grid -> cos^5(xyz) */
double *cos5xyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = pow(cos(x_(i)*y_(i)*z_(i)),5);
      
    in = fi;
  }
  return f;
}

/* f: grid -> sin(xyz) */
double *sinxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sin(x_(i)*y_(i)*z_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> sin^3(xyz) */
double *sin3xyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = pow(sin(x_(i)*y_(i)*z_(i)),3);
      
    in = fi;
  }
  return f;
}

/* f: grid -> cosh(xyz) */
double *coshxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = cosh(x_(i)*y_(i)*z_(i));
      
    in = fi;
  }
  return f;
}


/* f: grid -> sinh(xyz) */
double *sinhxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = sinh(x_(i)*y_(i)*z_(i));
      
    in = fi;
  }
  return f;
}


/* f: grid -> tanh(xyz) */
double *tanhxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = tanh(x_(i)*y_(i)*z_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> log(1+(xyz)^2) */
double *logxyz_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = log(1+pow(x_(i)*y_(i)*z_(i),2));
      
    in = fi;
  }
  return f;
}


/* f: grid -> x^2*cos(x*y^3*z) +y^3*z^4*x*sin(x^2*z) */
double *mix1_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = pow(x_(i),2)*cos(x_(i)*pow(y_(i),3)*z_(i))+
             pow(y_(i),3)*pow(z_(i),4)*x_(i)*sin(pow(x_(i),2)*z_(i));
      
    in = fi;
  }
  return f;
}

/* f: grid -> log(1+z^2)*cosh{sqrt(x^2+y^2+z^2)+sin(exp(x*y*z))} */
double *mix2_f(const Grid_T *const grid)
{
  double *f = alloc_double(grid->nn);
  unsigned in/* initial point */,fi/* final point */;
  unsigned i,pa;
  
  in = 0;
  FOR_ALL(pa,grid->patch)
  {
    Patch_T *patch = grid->patch[pa];
    fi = total_nodes_patch(patch) + in;
    
    for (i = in; i < fi; ++i)
      f[i] = log(1+pow(z_(i),2))*cosh(
        sqrt(pow(x_(i),2)+pow(y_(i),2)+pow(z_(i),2))+
          sin(exp(x_(i)*y_(i)*z_(i)))
                                );
      
    in = fi;
  }
  return f;
}



