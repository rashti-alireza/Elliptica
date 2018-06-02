/*
// Alireza Rashti
// May 2018
*/

#include "reading_input_file.h"

/* reading parameters */
int read_input_file(char *const path)
{
  FILE *input;// input file
  char *buff;
  
  null_pathEr(path);
  
  input = fopen(path,"r");
  pointerEr(input);
  
  /* parsing and reading input file and making buffer */
  buff = make_buffer(input);
  
  /* making parameters */
  populate_parameters(buff);
  
  /* printing parameters */
  if (check_print(PRINT_PARAMETER))
    print_parameter();
  
  /* clean up */
  fclose(input);
  free(buff);
  
  return EXIT_SUCCESS;
}

/* populating parameters:
// breaking down the buffer to a string like "par_l = par_r" and then
// add par_l and par_r to the parameter data base
*/
static void populate_parameters(char *const buff)
{
  char *tok, delimit[2] = {ENTER,'\0'}, delimit2[2] = {EQUAL,'\0'};
  char *buff2;
  
  buff2 = strdup(buff);
  
  tok = strtok(buff2,delimit);
  while (tok != 0)
  {
    char *tok2;
    char *par_l = 0, *par_r = 0;
    enum FLOW f = LEFT;
    
    tok2 = strtok(tok,delimit2);
    while (tok2 != 0)
    {
      
      if (f == LEFT)
      {
        par_l = strdup(tok2);
      }
      else if (f == RIGHT)
      {
        par_r = strdup(tok2);
      }
      
      f = RIGHT;
      tok2 = strtok(0,delimit2);
      
    }// while (tok2 != 0)
    
    add_parameter(par_l,par_r);
    
    tok = strtok(0,delimit);
    
    if (par_l != 0)	free(par_l);
    if (par_r != 0)	free(par_r);
    
  }// while (tok != 0)
  
}
  

/* parsing and reading input file and making buffer */
static void *make_buffer(FILE *input)
{
  char *buff = 0, c;
  int i;
  
  i = 0;
  while (!feof(input))
  {
    c = fgetc(input);
    
    /* excluding out the spaces */
    if (c == SPACE)
      continue;
    
    /* excluding out the comments */
    else if (c == COMMENT)
    {
      while (c != ENTER || c != EOF)
        c = fgetc(input);
    }
    
    else
    {
      buff = realloc(buff,i+1);
      pointerEr(buff);
      buff[i++] = c;
    }
    
  }// end of while

    buff = realloc(buff,i+1);
    pointerEr(buff);
    buff[i++] = END;
    
  return buff;
}
