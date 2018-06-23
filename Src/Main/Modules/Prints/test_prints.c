/*
// Alireza Rashti
// June 2018
*/

#include "test_prints.h"

/* printing different quantities for test */

/* print parameters */
void pr_parameters(void)
{
  FILE *f;
  char dir[1000]={'\0'}, *path;
  extern Parameter_T **parameters_global;
  int i = 0;
  Flag_T flg;
  
  path = get_parameter_value_S("output_dir",&flg);
  parameterEr(flg);
  sprintf(dir,"%s/parameters.out",path);
  f = fopen(dir,"w");
  pointerEr(f);
  
  fprintf(f,SECTION"Parameters"SECTION"\n");
  
  while (parameters_global[i] != 0)
  {
    fprintf(f,"%s = %s\n",parameters_global[i]->lv,parameters_global[i]->rv);
    i++;
  }
  
  fclose(f);
}

/* if print option is on return 1 otherwise 0 */
int test_print(Print_T f)
{
  char *on; 
  
  switch(f)
  {
    case PRINT_PARAMETER:
      on = get_parameter_value_S("print_parameter",0);
      if (on == 0) return 0;
      
      if (strcmp(on,"yes") == 0 || strcmp(on,"y") == 0)
        return 1;
        
      break;
      
    default:
      break;
  }
  
  return 0;
}
