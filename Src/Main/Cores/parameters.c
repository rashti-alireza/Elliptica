/*
// Alireza Rashti
// May 2018
*/

#include "parameters.h"

/* adding left value and right value to parameter data base*/
void add_parameter(char *lv, char *rv)
{
  assert(lv != 0);
  
  Parameter_T *par;
  int l = strlen(lv),
      r = strlen(rv);
  
  par = alloc_parameter(global_parameter);
  
  par->lv = malloc(l+1);
  checkup(par->lv);
  strcpy(par->lv,lv);
  
  par->rv = malloc(r+1);
  checkup(par->rv);
  strcpy(par->rv,rv);
}

/* having parameter name, it returns a pointer to 
// the corresponding parameter
*/
void *get_parameter(char *const par_name)
{
  int i;
  
  i = 0;
  while (global_parameter[i] != 0)
  {
    if (strcmp(global_parameter[i]->lv,par_name) == 0)
      return global_parameter[i];
  }
  
  return 0;
}

/* having the parameter name and its kind - NUMERIC or LITERAL-,
// it returns the value of parameter.
// note: if it is NUMERIC the value is written in value otherwise
// the value is returned.
*/
void *get_parameter_value(char *const par_name,Flag_T kind, double *value)
{
  extern Parameter_T **global_parameter;
  char *v
  int i;
  
  i = 0;
  while (global_parameter[i] != 0)
  {
    if (strcmp(global_parameter[i]->lv,par_name) == 0)
    {
      if (kind == NUMERIC)
      {
        *value = strtod(global_parameter[i]->rv,0);
        return 0;
      }
      else if (kind == LITERAL)
      {
        return global_parameter[i]->rv;
      }
      else
        bad_input();
    }
    
    i++;
  }
  
  return 0;
}
