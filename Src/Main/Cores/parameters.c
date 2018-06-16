/*
// Alireza Rashti
// May 2018
*/

#include "parameters.h"

/* adding left value and right value to parameter data base*/
void add_parameter(char *lv, char *rv)
{
  pointerEr(lv);
  
  Parameter_T *par;
  
  par = get_parameter(lv);
  if (par)
    abortEr_s("This parameter \"%s\" has already been added!\n",lv);
    
  par = alloc_parameter(&parameters_global);
  
  par->lv = strdup(lv);
  par->rv = strdup(rv);
}

/* having parameter name, it returns a pointer to 
// the corresponding parameter
*/
void *get_parameter(char *const par_name)
{
  int i;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp(parameters_global[i]->lv,par_name) == 0)
      return parameters_global[i];
    i++;
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
  int i;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp(parameters_global[i]->lv,par_name) == 0)
    {
      if (kind == NUMERIC)
      {
        *value = strtod(parameters_global[i]->rv,0);
        return 0;
      }
      else if (kind == LITERAL)
      {
        return parameters_global[i]->rv;
      }
      else
        bad_inputEr();
    }
    
    i++;
  }
  
  return 0;
}
