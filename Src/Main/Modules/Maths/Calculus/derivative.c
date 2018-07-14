/*
// Alireza Rashti
// July 2018
*/

#include "derivative.h"
#define DELIMIT '|'

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
double *Df(const Field_T *const f,const char *task)
{
  if (!task)
    abortEr("The task is empty!\n");
    
  double *r = 0;
  const char *der_par = get_parameter_value_S("Derivative_Method",0);
  Derivative_T Der_e = derivative_type(der_par,task);
  
  return r;
}

/* getting both task and par, it finds the derivative method.
// if task has a derivative method, it is preferred over par.
// ->return value: derivative method in enum format.
*/
static Derivative_T derivative_type(const char *const par,const char *const task)
{
  const char delimit = DELIMIT;
  Derivative_T type = UNDEFINED_DERIVATIVE;
  char *s = dup_s(task);
  char *rs = 0;
  
  tok_s(s,delimit,&rs);
  type = str2enum(rs);
  if (s) free(s);
  
  /* if check parameter if no info is in task */
  if (type == UNDEFINED_DERIVATIVE)
  {
    type = str2enum(par);
  } 
  
  /* if still no type has been found => no info in parameter and task */
  if (type == UNDEFINED_DERIVATIVE)
    abortEr("No Derivative Method is defined "
      "in parameter file or in function task.\n");
  
  return type;
}

/* getting a derivative method in string format and 
// returing the Derivative_T.
// ->return value: Derivative_T and if not found UNDEFINED_DERIVATIVE.
*/
static Derivative_T str2enum(const char *const str)
{
  Derivative_T type = UNDEFINED_DERIVATIVE;
  
  if (strcmp_i(str,"Spectral"))
    type = SPECTRAL;
  else if (strcmp_i(str,"Finite_Difference"))
    type = FINITE_DIFF;
  
  return type;
}
