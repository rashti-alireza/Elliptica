/*
// Alireza Rashti
// June 2018
*/

#include "test_prints.h"

/* printing different quantities for test */

/* print parameters */
void print_parameters(void)
{
  FILE *f;
  char dir[1000]={'\0'};
  extern Parameter_T **global_parameter;
  int i = 0;
  
  sprintf(dir,"%s/parameters.out",global_path);
  f = fopen(dir,"w");
  pointerEr(f);
  
  while (global_parameter[i] != 0)
  {
    fprintf(f,"%s = %s \n",global_parameter[i]->lv,global_parameter[i]->rv);
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
      on = get_parameter_value("print_parameter",LITERAL,0);
      
      if (strcmp(on,"yes") || strcmp(on,"1"))
        return 1;
        
      break;
      
    default:
      break;
  }
  
  return 0;
}