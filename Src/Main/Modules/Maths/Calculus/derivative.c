/*
// Alireza Rashti
// July 2018
*/

#include "derivative.h"
#define DELIMIT '|'
#define COMMA ','

/* derivative function. it does the task and return the result.
// note: this function allocate memory for the result.
//
// list of available derivative methods:
// =====================================
// o. spectral = for spectral method
//
// task format example:
// ===================
//
// Cartesian x derivative => *task = "x"
// Cartesian z followed by Cartesian y derivative => *task = "z,y"
// Cartesian xx derivative => *task = "x,x"
// and so forth. furthermore for derivative in curvilinear coords we have:
// Curvilinear a derivative followed by b followed by c => *task ="a,b,c"
// if one wants to override the default derivative method defined in 
// the input file, they can append the task by "|derivative type"; e.g.
// Cartesian x derivative with finite difference:
// => *task = "x DELIMIT Finite_Difference", DELIMIT is a macro defined above.
// so for example if DELIMIT is | then *task = "x | Finite_Difference".
// ->return value: derivative of Field_T f accordingly, null for error.
*/
double *Df(Field_T *const f,const char *task)
{
  /* check up */
  if (!f)
    abortEr("The field is empty!\n");
  if (!f->values && !f->coeffs)  
    abortEr("The field is empty!\n");
  if (!task)
    abortEr("The task is empty!\n");

  double *r = 0;
  const char *der_par = get_parameter_value_S("Derivative_Method",0);
  unsigned Ndir;
  Method_T method_e = derivative_method(der_par,task);
  Direction_T  *dir_e = derivative_direction(task,&Ndir);
  
  if (method_e == SPECTRAL)
    r = take_spectral_derivative(f,dir_e,Ndir);
  else
    abortEr("There is no such derivative method defined for this function.\n");
    
  free(dir_e);
  return r;
}

/* taking spectral derivative of the given field 
// based on direction of derivative.
// ->return value: values of the resultant.
*/
static double *take_spectral_derivative(Field_T *const f,const Direction_T  *const dir_e,const unsigned Ndir)
{
  double *deriv = 0;
  Field_T *ff[2];
  unsigned bck,frd;
  unsigned i;
  
  /* 3-D */
  if (f->dim == 3)
  {
    
    ff[0] = init_field_3d("tmp1",f->grid);
    ff[1] = init_field_3d("tmp2",f->grid);
    
    ff[0]->values = spectral_derivative_1d(f,dir_e[i]);
    
    for (i = 1; i < Ndir; ++i)
    {
      frd = i%2;
      bck = i/2;
      
      ff[frd]->values = spectral_derivative_1d(ff[bck],dir_e[i]);
      /* make next ff ready */
      free(ff[bck]->values);
      free(ff[bck]->coeffs);
      ff[bck]->values = 0;
      ff[bck]->coeffs = 0;
      //THINK ABOUT THIS!!!
    }
    
    /* free leftovers */
    free_field(ff[0]);
    free_field(ff[1]);
  }/* end of if (f->dim == 3) */
  else
    abortEr("No such Dimension is defined for this function.\n");
  
  return deriv;
}

/* finding all of types of derivatives and put them into 
// an array of Direction_T with size n in order that they are written.
// note: this function allocate memory.
// ->return value: array of Direction_T and number of this arrays
*/
static Direction_T *derivative_direction(const char *const task,unsigned *const n)
{
  Direction_T *e = 0;
  char *savestr,*str = dup_s(task);
  char *tok = tok_s(str,DELIMIT,&savestr);
  
  if (!tok)
    abortEr_s("There is No direction in %s.\n",task);
  
  *n = 0;  
  while (tok)
  {
    e = realloc(e,(*n+1)*sizeof(*e));
    pointerEr(e);
    e[*n] = str2enum_direction(tok);
    tok = tok_s(0,COMMA,&savestr);
    ++(*n);
  }
  
  free(str);
  
  return e;
}

/* getting both task and par, it finds the derivative method.
// if task has a derivative method, it is preferred over par.
// ->return value: derivative method in enum format.
*/
static Method_T derivative_method(const char *const par,const char *const task)
{
  const char delimit = DELIMIT;
  Method_T type = UNDEFINED_METHOD;
  char *s = dup_s(task);
  char *rs = 0;
  
  tok_s(s,delimit,&rs);
  type = str2enum_method(rs);
  if (s) free(s);
  
  /* if check parameter if no info is in task */
  if (type == UNDEFINED_METHOD)
  {
    type = str2enum_method(par);
  } 
  
  /* if still no type has been found => no info in parameter and task */
  if (type == UNDEFINED_METHOD)
    abortEr("No Derivative Method is defined "
      "in parameter file or in function task.\n");
  
  return type;
}

/* getting a derivative method in string format and 
// returing the Method_T.
// ->return value: Method_T and if not found UNDEFINED_Method.
*/
static Method_T str2enum_method(const char *const str)
{
  Method_T type = UNDEFINED_METHOD;
  
  if (strcmp_i(str,"Spectral"))
    type = SPECTRAL;
  else if (strcmp_i(str,"Finite_Difference"))
    type = FINITE_DIFF;
  
  return type;
}

/* getting a derivative direction in string format and 
// returning the Direction_T.
// ->return value: Direction_T and if not found error.
*/
static Direction_T str2enum_direction(const char *const str)
{
  if (strcmp_i(str,"x"))
    return x_DIR;
  else if (strcmp_i(str,"y"))
    return y_DIR;
  else if (strcmp_i(str,"z"))
    return z_DIR;
  else if (strcmp_i(str,"a"))
    return a_DIR;
  else if (strcmp_i(str,"b"))
    return b_DIR;
  else if (strcmp_i(str,"c"))
    return c_DIR;
  else
    abortEr_s("There is no derivative such as %s.\n",str);
  
  return UNDEFINED_DIR;
}

/* taking 3-D spectral derivative in specified direction
// ->return value: derivative.
*/
static double *spectral_derivative_1d(Field_T *const f,const Direction_T dir_e)
{
  double *der = 0;
  double *df_dN[3];
  Grid_T *const grid = f->grid;
  const Patch_T *patch;
  unsigned pa,d;
  
  FOR_ALL(pa,grid->patch)
  {
    patch = grid->patch[pa];
    
    for (d = 0; d < 3; ++d)
    {
      if(patch->basis[d]       == Chebyshev_Tn_BASIS && 
         patch->collocation[d] == Chebyshev_Extrema)
        df_dN[d] = derivative_Chebyshev_Tn_1d(f,patch,d);
      else
        abortEr("There is no such basis or collocation defined for this function.\n");
    }
    
    //JACOBIAN MULTIPLICATION
  }
  
  return der;
}

/* taking derivative for f(x)= sum_{0}^{n-1}a_n*T_n(x), in the specified
// patch.
// note: it is 1-D, and x is in [-1,1]. for 3-d one needs to combine these
// with Jacobian transformation.
// ->return value: df(x)/dx.
*/
static double *derivative_Chebyshev_Tn_1d(Field_T *const f,const Patch_T *const patch,const unsigned dir)
{
  const double *const coeffs = make_coeffs_1d(f,dir,patch);
  const unsigned *const n = patch->n;
  const unsigned nc = patch->nc;
  const unsigned nn = total_nodes_patch(patch);
  const unsigned B = n[dir]-1;
  double *der = alloc_double(nn);
  double *x = make_normalized_collocation_1d(patch,dir);
  unsigned l;
  
  #pragma omp parallel for
  for (l = 0; l < nn; ++l)
  {
    unsigned m;
    unsigned L = nc+l;
    for (m = 2; m < B; ++m)
    {
      unsigned lc = L_c(l,m,n,dir);
      double u = Cheb_Un(m-1,x[m]);
      der[L] += m*coeffs[][lc]*u;
    }
    der[L] += coeffs[][L_c(l,1,n,dir)];
    der[L] *= 2;
  }
  
  free(x);
}
