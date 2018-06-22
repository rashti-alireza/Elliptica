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
Parameter_T *get_parameter(char *const par_name)
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
char *get_parameter_value(char *const par_name,Flag_T kind, double *value)
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

/* reading the input file and make all of parameters 
// with their defualt value to the data base 
*/
int make_parameters(char *const path)
{
  char folder[100]={'\0'}, *name;
  char *path2;
  
  read_input_file(path);
  
  /* setting the default value of parameters if they are needed and
  not provided by the inputfile */
  set_default_parameter();
  
  /* making a folder at the directory of 
  // input file with the name of "inputfile_output"
  // or with the given name in input file 
  // and rewritting global_path with new directory path 
  */
  name = get_parameter_value("output_directory_name",LITERAL,0);
  sprintf(folder,"%s_output",name);
  path2 = make_directory(path_global,folder,YES);
  free(path_global);
  path_global = path2;
  
  add_parameter("output_path",path_global);
  
  /* printing parameters */
  if (test_print(PRINT_PARAMETER))
    pr_parameters();
  
  return EXIT_SUCCESS;
}

