/*
// Alireza Rashti
// May 2018
*/

#include "parameters.h"

/* adding left value and right value to parameter data base*/
void add_parameter(const char *const lv, const char *const rv)
{
  pointerEr(lv);
  
  Parameter_T *par;
  
  par = get_parameter(lv);
  if (par)
    abortEr_s("This parameter \"%s\" has already been added!\n",lv);
    
  par = alloc_parameter(&parameters_global);
  par->lv = dup_s(lv);
  if (!strcmp(rv,"") || !strcmp(rv," ")) par->rv = 0;
  else par->rv = dup_s(rv);
}

/* having parameter name, it returns a pointer to 
// the corresponding parameter
*/
Parameter_T *get_parameter(const char *const par_name)
{
  int i;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
      return parameters_global[i];
    i++;
  }
  
  return 0;
}

/* having the parameter name,
// it returns the INTEGER value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: integer value of parameter.
*/
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  int v = INT_MAX;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = atoi(parameters_global[i]->rv);
      f = FOUND;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter %s couldn't be found.\n",par_name,file,line);
  }

  return v;
}

/* having the parameter name,
// it returns the DOUBLE value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: double value of parameter.
*/
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  double v = DBL_MAX;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = strtod(parameters_global[i]->rv,0);
      f = FOUND;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter %s couldn't be found.\n",par_name,file,line);
  }

  return v;
}

/* having the parameter name,
// it returns the STRING value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: string value of parameter.
*/
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  char *v = 0;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = parameters_global[i]->rv;
      f = FOUND;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter %s couldn't be found.\n",par_name,file,line);
  }
      
  return v;
}

/* reading the input file and make all of parameters 
// with their default value to the data base 
*/
int make_parameters(const char *const path)
{
  char folder[100]={'\0'};
  const char *name;
  char *new_path;
  
  read_input_file(path);
  
  /* setting the default value of parameters if they are needed and
  not provided by the inputfile */
  set_default_parameter();
  
  /* making a folder at the directory of 
  // input file with the name of "inputfile_output"
  // or with the given name in input file 
  // and rewriting global_path with new directory path 
  */
  name = GetParameterS_E("output_directory_name");
  sprintf(folder,"%s_output",name);
  new_path = make_directory(path_global,folder);
  free(path_global);
  path_global = dup_s(new_path);
  add_parameter("output_directory_path",new_path);
  
  /* printing parameters */
  if (test_print(PRINT_PARAMETERS))
    pr_parameters();
  
  free(new_path);
  
  return EXIT_SUCCESS;
}
