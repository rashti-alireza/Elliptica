/*
// Alireza Rashti
// May 2018
*/

#include "reading_input_file.h"

/* reading input file and make appropriate directories for output files */
int read_input_file(char *const path)
{
  FILE *input;// input file
  char *buff;
  char name[100]={'\0'};
  char *path2;
  
  null_pathEr(path);
  
  input = fopen(path,"r");
  null_pathEr(input);  
  
  /* making a folder at the directory of 
  // input file with the name of "inputfile_output" 
  // and rewritting global_path with new directory path 
  */
  sprintf(name,"%s_output",inputfile_name_global);
  path2 = make_directory(path_global,name,YES);
  free(path_global);
  path_global = path2;
  
  add_parameter("output_dir",path_global);
  
  /* parsing and reading input file and making buffer */
  buff = make_buffer(input);
  
  //TEST_START
    //printf("buff:\n%s\n",buff);
  //end
  
  /* making parameters */
  populate_parameters(buff);
  
  /* printing parameters */
  if (test_print(PRINT_PARAMETER))
    print_parameters();
  
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
  char *tok,*savestr;
  char delimit[2] = {ENTER,'\0'}, delimit2[2] = {EQUAL,'\0'};
  char *buff2;
  
  buff2 = strdup(buff);
  
  tok = strtok_r(buff2,delimit,&savestr);
  while (tok != 0)
  {
    char *subtok;
    char *par_l = 0, *par_r = 0;
    enum FLOW f = LEFT;
    
    subtok = strtok(tok,delimit2);
    while (subtok != 0)
    {
      
      if (f == LEFT)
      {
        par_l = strdup(subtok);
      }
      else if (f == RIGHT)
      {
        par_r = strdup(subtok);
      }
      
      f = RIGHT;
      subtok = strtok(0,delimit2);
      
    }// while (subtok != 0)
    
    if (par_l != 0 && par_r != 0)
      add_parameter(par_l,par_r);
    
    tok = strtok_r(0,delimit,&savestr);
  
    //TEST_START
      //printf("par_l = %s, par_r = %s, tok = %s\n",par_l,par_r,tok);
    //end
    
    if (par_l != 0)	free(par_l);
    if (par_r != 0)	free(par_r);
    
  }// while (tok != 0)
  
  free(buff2);
}
  

/* parsing and reading input file and making buffer */
static void *make_buffer(FILE *input)
{
  char *buff = 0, c;
  int i;
  long j;
  
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
        //test
        //printf("%c",c);
        //end
      }
      
      if (c == ENTER)
      {
        buff = realloc(buff,i+1);
        pointerEr(buff);
        buff[i++] = c;
      }
    }
    
    else
    {
      buff = realloc(buff,i+1);
      pointerEr(buff);
      buff[i++] = c;
    }
    
    if (j > MAX_NUM_CHAR)
      abortEr("Input file has been written wrong!\n"
        "It might beacuse either it is too long or "
          "beacuse its format has not been met\n");
    j++;      
  }// end of while

    buff = realloc(buff,i+1);
    pointerEr(buff);
    buff[i] = END;
    
  return buff;
}
