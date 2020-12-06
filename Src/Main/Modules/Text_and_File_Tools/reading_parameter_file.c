/*
// Alireza Rashti
// May 2018
*/

#include "reading_parameter_file.h"

/* reading input file and make appropriate directories for output files */
void read_input_file(const char *const path)
{
  FILE *input;/* input file */
  char *buff;
  
  null_pathEr(path);
  
  input = Fopen(path,"r");
  null_pathEr(input);  
  
  /* parsing and reading input file and making buffer */
  buff = make_buffer(input);
  
  /* is file empty? */
  if (strlen(buff) < 2)
    Error0("The parameter file is empty.\n");
  
  //TEST_START
    //printf("buff:\n%s\n",buff);
  //end
  
  /* populating parameters */
  populate_parameters(buff);
  
  /* clean up */
  Fclose(input);
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
    enum FLOW f = e_Left;
    
    subtok = tok_s(tok,delimit2,&savestr2);
    while (subtok != 0)
    {
      
      if (f == e_Left)
      {
        par_l = dup_s(subtok);
      }
      else if (f == e_Right)
      {
        par_r = dup_s(subtok);
      }
      
      f = e_Right;
      subtok = tok_s(0,delimit,&savestr2);
      
    }/* while (subtok != 0) */
    
    if (par_l)
    {
      if (!par_r)
        Errors("No value for parameter '%s'\n",par_l);
      add_parameter(par_l,par_r);
    }
    
    tok = tok_s(0,delimit,&savestr);
  
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
  Uint i;
  Uint long j;
  
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
        IsNull(buff);
        buff[i++] = (char)c;
      }
    }
    
    /* joining broken lines */
    else if (c == BACK_SLASH)
    {
      c = fgetc(input);
      if (c != ENTER)/* if after \ is not enter */
        Error0("No enter after '\\'. Broken line is incorrect.\n");
        
      c = fgetc(input);/* skip enter */
      while (c != ENTER && c != EOF)
      {
        if (c == SPACE ||  c == TAB)/* trim white spaces */
        {
          c = fgetc(input);
          continue;
        }
        else if (c == BACK_SLASH)/* if multi broken line */
        {
          c = fgetc(input);
          if (c != ENTER)/* if after \ is not enter */
            Error0("No enter after '\\'. Broken line is incorrect.\n");
          c = fgetc(input);/* skip enter */
          continue;
        }
        //else if (c == EQUAL)
          //Error0("It found '=' sign in broken line. Broken line is incorrect.\n");
        else if (c == COMMENT)/* time off comment */
        {
          while (c != ENTER && c != EOF) c = fgetc(input);
          if (c == ENTER)                break;
        }
        
        buff = realloc(buff,(i+1)*sizeof(*buff));
        IsNull(buff);
        buff[i++] = (char)c;
        c = fgetc(input);
      }/* end of while (c != ENTER && c != EOF) */
      if (c == ENTER)
      {
        buff = realloc(buff,(i+1)*sizeof(*buff));
        IsNull(buff);
        buff[i++] = (char)c;
      }
    }/* end of else if (c == BACK_SLASH) */
    else
    {
      buff = realloc(buff,(i+1)*sizeof(*buff));
      IsNull(buff);
      buff[i++] = (char)c;
    }
    
    if (j > MAX_NUM_CHAR)
      Error0("Input file has been written wrong!\n"
        "It might beacuse either it is too long or "
          "beacuse its format has not been met\n");
    j++;      
  }/* end of while */

  //buff = realloc(buff,(i+1)*sizeof(*buff));
  //IsNull(buff);
  //buff[i] = END;
  buff[i-1] = END;
  return buff;
}
