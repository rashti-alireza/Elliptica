/*
// Alireza Rashti
// May 2018
*/

#include "checkup.h"
#define MAX_ARR 600
#define MAX_WIDTH 73

/* checking up if the allocation of the pointer is failed;
// if it is failed then exit.
*/
void checkup_pointer_error(const void *const p, const char *const file, const int line)
{
  if (p == NULL)
  {
    pr_line_stderr('x');
    fprintf(stderr,ERROR_MASSAGE
      "\"Pointer allocation was failed.\"\n"
        "File: %s\nLine: %d\n",file,line);
    pr_line_stderr('x');
    abort();
  }
}

/* bad input */
void bad_input_error(const char *const file, const int line)
{
  pr_line_stderr('x');
  fprintf(stderr,ERROR_MASSAGE
          "\"There is no such input for the function.\"\n"
          "File: %s\nLine: %d\n",file,line);
  pr_line_stderr('x');
  abort();
}

/* null path directory */
void null_path_error(const void *const path,const char *const file, const int line)
{
  if (path == NULL)
  {
    pr_line_stderr('x');
    fprintf(stderr,ERROR_MASSAGE
      "\"The directory path is null.\"\n"
        "File: %s\nLine: %d\n",file,line);
    pr_line_stderr('x');
    abort();
  }
}

/* general purpose error, if pr_strdout == 1 it ALSO prints in standard output */
void abort_error(const char *const massage,const char *const file, const int line,const int pr_stdout)
{
  if (pr_stdout)
  {
    printf(ERROR_MASSAGE"%s\nFile: %s\nLine: %d\n\n",massage,file,line);
    fflush(stdout);
  }
  pr_line_stderr('x');
  fprintf(stderr,ERROR_MASSAGE"%s\nFile: %s\nLine: %d\n",massage,file,line);
  pr_line_stderr('x');
  abort();
}

/* general purpose error with more detail string version */
void abort_error_string(const char *const massage,const char *const detail,const char *const file, const int line)
{
  char msg[MAX_ARR] = {'\0'};
  
  sprintf(msg,massage,detail);
  
  pr_line_stderr('x');
  fprintf(stderr,ERROR_MASSAGE"%s\nFile: %s\nLine: %d\n",msg,file,line);
  pr_line_stderr('x');
  abort();
}

/*chech if the parameter has been found */
void check_parameter(const Flag_T flg,const char *const file, const int line)
{
  if (flg != FOUND)
  {
    pr_line_stderr('x');
    fprintf(stderr,ERROR_MASSAGE
      "\"Parameter was not found.\"\n"
        "File: %s\nLine: %d\n",file,line);
    pr_line_stderr('x');
    abort();
  }
  
}

/* print a line in stderr with char c */
static void pr_line_stderr(const char c)
{
  int i;
  
  fflush(stdout);
  for (i = 0; i < MAX_WIDTH; i++)
    fprintf(stderr,"%c",c);
  fprintf(stderr,"\n");
  
  fflush(stderr);
}
