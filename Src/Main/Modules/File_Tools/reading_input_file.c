/*
// Alireza Rashti
// May 2018
*/

#include "reading_input_file.h"

/* reading input file and make appropriate directories for output files */
void read_input_file(const char *const path)
{
  FILE *input;/* input file */
  char *buff;
  
  null_pathEr(path);
  
  input = fopen(path,"r");
  null_pathEr(input);  
  
  /* parsing and reading input file and making buffer */
  buff = make_buffer(input);
  
  /*TEST_START
    //printf("buff:\n%s\n",buff);
  //end */
  
  /* populating parameters */
  populate_parameters(buff);
  
  /* clean up */
  fclose(input);
  free(buff);

}

/* populating parameters:
// breaking down the buffer to a string like "par_l = par_r" and then
// add par_l and par_r to the parameter data base
*/
static void populate_parameters(const char *const buff)
{
  char *tok,*savestr,*savestr2;
  char delimit = ENTER, delimit2 = EQUAL;
  char *buff2;
  
  buff2 = dup_s(buff);
  
  tok = tok_s(buff2,delimit,&savestr);
  while (tok != 0)
  {
    char *subtok;
    char *par_l = 0, *par_r = 0;
    enum FLOW f = LEFT;
    
    subtok = tok_s(tok,delimit2,&savestr2);
    while (subtok != 0)
    {
      
      if (f == LEFT)
      {
        par_l = dup_s(subtok);
      }
      else if (f == RIGHT)
      {
        par_r = dup_s(subtok);
      }
      
      f = RIGHT;
      subtok = tok_s(0,delimit2,&savestr2);
      
    }/* while (subtok != 0) */
    
    if (par_l != 0 && par_r != 0)
      add_parameter(par_l,par_r);
    
    tok = tok_s(0,delimit,&savestr);
  
    /*TEST_START
      //printf("par_l = %s, par_r = %s, tok = %s\n",par_l,par_r,tok);
    //end */
    
    if (par_l != 0)	free(par_l);
    if (par_r != 0)	free(par_r);
    
  }/* while (tok != 0) */
  
  free(buff2);
}
  

/* parsing and reading input file and making buffer */
static void *make_buffer(FILE *const input)
{
  char *buff = 0;
  int c;
  unsigned i;
  unsigned long j;
  
  j = 0;
  i = 0;
  while (!feof(input))
  {
    c = fgetc(input);
    
    /* excluding out the spaces and tabs*/
    if (c == SPACE ||  c == TAB)
      continue;
    
    /* excluding out the comments */
    else if (c == COMMENT)
    {
      while (c != ENTER && c != EOF)
      {
        c = fgetc(input);
      }
      
      if (c == ENTER)
      {
        buff = realloc(buff,(i+1)*sizeof(*buff));
        pointerEr(buff);
        buff[i++] = (char)c;
      }
    }
    
    else
    {
      buff = realloc(buff,(i+1)*sizeof(*buff));
      pointerEr(buff);
      buff[i++] = (char)c;
    }
    
    if (j > MAX_NUM_CHAR)
      abortEr("Input file has been written wrong!\n"
        "It might beacuse either it is too long or "
          "beacuse its format has not been met\n");
    j++;      
  }/* end of while */

    buff = realloc(buff,(i+1)*sizeof(*buff));
    pointerEr(buff);
    buff[i] = END;
    
  return buff;
}
