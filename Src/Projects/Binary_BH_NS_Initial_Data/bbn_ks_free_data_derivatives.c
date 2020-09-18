#include "bbn_ks_free_date_analytic.h"

#define MAX_ARR (100)

static void 
function_name
  (
  char *const fname/* function name */,
  const char *const stem/* eg. bbn_ks_X */,
  const char *const derive/* x,y,z */
  );

/* interpret analytic derivative */
double bbn_ks_derivative KS_deriv_func_args_macro
{
  unsigned f;
  char fname[MAX_ARR] = {'\0'};
  
  /* realize the function name */
  function_name(fname,stem,derivs);
  
  /* find the index in the data base and call */
  while (derive_name_db[f])
  {
    if (!strcmp(fname,derive_name_db[f]))
      return derive_func_db[f]KS_func_pass_args_macro;
  }
  
  /* it must not reach here! */
  Error0("No such function found.\n");
  
  return DBL_MAX;
}

/* given stem and derives it finds the correspondence derivative
// function name and writes it in fname. */ 
static void 
function_name
  (
  char *const fname/* function name */,
  const char *const stem/* eg. bbn_ks_X */,
  const char *const derive/* x,y,z */
  )
{
  char func[MAX_ARR] = {'\0'};
  sprintf(fname,"%s",func);
}

