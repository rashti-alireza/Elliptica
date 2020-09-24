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
  f = 0;
  while (derive_name_db[f])
  {
    if (!strcmp(fname,derive_name_db[f]))
      return derive_func_db[f]KS_func_pass_args_macro;
    f++;
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
  const char *const stem/* eg. bbn_ks_X(x,y,z) */,
  const char *const derive/* x,y,z */
  )
{
  char der[MAX_ARR] = {'\0'};
  char root[MAX_ARR] = {'\0'};
  char d[MAX_ARR] = {'\0'};
  char s[MAX_ARR] = {'\0'};
  char *ps;
  unsigned i,j;
  
  /* trim (x,y,z) */
  sprintf(s,"%s",stem);
  ps = strchr(s,'(');
  ps[0] = '\0';
  
  /* trim prefix "bbn_ks_" */
  ps = s+strlen(KS_prefix);
  sprintf(root,"%s",ps);
  
  /* remove spaces */
  i = j = 0;
  while(derive[i] != '\0')
  {
    if (derive[i] != ' ')
    {
      s[j] = derive[i];
      j++;
    }
    i++;
  }
  s[j] = '\0';
  
  /* 3rd order (x,y,z) */
  if (strstr(s,"x") && 
      strstr(s,"y") &&
      strstr(s,"z"))
  {
    sprintf(der,"D0D1D2");
    sprintf(d,"ddd");
  }
  /* 3rd order:
  //  (x,(z,2))|(x,(y,2))|(y,(x,2))|(y,(z,2))|(z,(x,2))|(z,(y,2)) */
  else if (!strcmp(s,"(x,(z,2))"))
  {
    sprintf(der,"D0D2D2");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(x,(y,2))"))
  {
    sprintf(der,"D0D1D1");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(y,(x,2))"))
  {
    /* symmetry */
    sprintf(der,"D0D1D0");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(y,(z,2))"))
  {
    sprintf(der,"D1D2D2");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(z,(x,2))"))
  {
    /* symmetry */
    sprintf(der,"D0D2D0");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(z,(y,2))"))
  {
    /* symmetry */
    sprintf(der,"D1D2D1");
    sprintf(d,"ddd");
  }
  /* (x,3) - (y,3) - (z,3) */
  else if (!strcmp(s,"(x,3)"))
  {
    /* symmetry */
    sprintf(der,"D0D0D0");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(y,3)"))
  {
    /* symmetry */
    sprintf(der,"D1D1D1");
    sprintf(d,"ddd");
  }
  else if (!strcmp(s,"(z,3)"))
  {
    /* symmetry */
    sprintf(der,"D2D2D2");
    sprintf(d,"ddd");
  }
  /* 2nd order (x,y) - (x,z) -(y,z) */
  else if (strstr(s,"x") && 
           strstr(s,"y"))
  {
    sprintf(der,"D0D1");
    sprintf(d,"dd");
  }
  else if (strstr(s,"x") && 
           strstr(s,"z"))
  {
    sprintf(der,"D0D2");
    sprintf(d,"dd");
  }
  else if (strstr(s,"y") && 
           strstr(s,"z"))
  {
    sprintf(der,"D1D2");
    sprintf(d,"dd");
  }
  /* 2rd order:
  //  (x,2)|(y,2)|(z,2) */
  else if (!strcmp(s,"(x,2)")) 
  {
    sprintf(der,"D0D0");
    sprintf(d,"dd");
  }
  else if (!strcmp(s,"(y,2)")) 
  {
    sprintf(der,"D1D1");
    sprintf(d,"dd");
  }
  else if (!strcmp(s,"(z,2)")) 
  {
    sprintf(der,"D2D2");
    sprintf(d,"dd");
  }
  /* 1st order x - y - z */
  else if (!strcmp(s,"x"))
  {
    sprintf(der,"D0");
    sprintf(d,"d");
  }
  else if (!strcmp(s,"y"))
  {
    sprintf(der,"D1");
    sprintf(d,"d");
  }
  else if (!strcmp(s,"z"))
  {
    sprintf(der,"D2");
    sprintf(d,"d");
  }
  else
    Error1("'%s' not recognized!",s);
  
  /* function name */
  sprintf(fname,"%s%s%s_%s",KS_prefix,d,root,der);
  /* test */
  // printf("'%s','%s' => %s\n",stem,derive,fname);
}




