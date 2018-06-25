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
  if (!strcmp(rv,"") || !strcmp(rv," ")) par->rv = 0;
  else par->rv = strdup(rv);
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
    if (strcmp_i(parameters_global[i]->lv,par_name))
      return parameters_global[i];
    i++;
  }
  
  return 0;
}

/* having the parameter name,
// it returns the INTEGER value of parameter.
// note:  if the parameter is found the flg gets FOUND
// value otherwise NONE value; unless flg = 0 which in this case
// it won't get any value
*/
int get_parameter_value_I(char *const par_name,Flag_T *flg)
{
  int v;
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
  
  if (flg != 0)
    *flg = f;
      
  return v;
}

/* having the parameter name,
// it returns the DOUBLE value of parameter.
// note:  if the parameter is found the flg gets FOUND
// value otherwise NONE value; unless flg = 0 which in this case
// it won't get any value
*/
double get_parameter_value_D(char *const par_name,Flag_T *flg)
{
  double v;
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
  
  if (flg != 0)
    *flg = f;
        
  return v;
}


/* having the parameter name,
// it returns the STRING value of parameter.
// note:  if the parameter is found the flg gets FOUND
// value otherwise NONE value; unless flg = 0 which in this case
// it won't get any value
*/
char *get_parameter_value_S(char *const par_name,Flag_T *flg)
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
  
  if (flg != 0)
    *flg = f;
      
  return v;
}

/* reading the input file and make all of parameters 
// with their defualt value to the data base 
*/
int make_parameters(char *const path)
{
  char folder[100]={'\0'}, *name;
  Flag_T flg;
  
  read_input_file(path);
  
  /* setting the default value of parameters if they are needed and
  not provided by the inputfile */
  set_default_parameter();
  
  /* making a folder at the directory of 
  // input file with the name of "inputfile_output"
  // or with the given name in input file 
  // and rewritting global_path with new directory path 
  */
  name = get_parameter_value_S("output_directory_name",&flg);
  parameterEr(flg);
  sprintf(folder,"%s_output",name);
  make_directory(&path_global,folder,YES);
  
  add_parameter("output_directory_path",path_global);
  
  /* printing parameters */
  if (test_print(PRINT_PARAMETER))
    pr_parameters();
  
  return EXIT_SUCCESS;
}
