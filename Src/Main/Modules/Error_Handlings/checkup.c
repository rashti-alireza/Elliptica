/*
// Alireza Rashti
// May 2018
*/

#include "checkup.h"
#define MAX_ARR 600

/* checking up if the allocation of the pointer is failed;
// if it is failed then exit.
*/
void checkup_pointer_error(const void *const p, const char *const file, const int line)
{
  if (p == NULL)
  {
    fprintf(stderr,ERROR_MASSAGE
      "\"Pointer allocation was failed.\"\n"
        "File: %s\nLine:%d\n",file,line);
    fflush(stderr);
    abort();
  }
}

/* bad input */
void bad_input_error(const char *const file, const int line)
{
    fprintf(stderr,ERROR_MASSAGE
      "\"There is no such input for the function.\"\n"
        "File: %s\nLine:%d\n",file,line);
    fflush(stderr);
    abort();
}

/* null path directory */
void null_path_error(const void *const path,const char *const file, const int line)
{
  if (path == NULL)
  {
    fprintf(stderr,ERROR_MASSAGE
      "\"The directory path is null.\"\n"
        "File: %s\nLine:%d\n",file,line);
    fflush(stderr);
    abort();
  }
}

/* general purpose error */
void abort_error(const char *const massage,const char *const file, const int line)
{
  
  fprintf(stderr,ERROR_MASSAGE"\"%s\""
        "File: %s\nLine:%d\n",massage,file,line);
  fflush(stderr);
  abort();
}

/* general purpose error with more detail string version */
void abort_error_string(const char *const massage,const char *const detail,const char *const file, const int line)
{
  char msg[MAX_ARR] = {'\0'};
  
  sprintf(msg,massage,detail);
  
  fprintf(stderr,ERROR_MASSAGE"%s"
      "File: %s\nLine:%d\n",msg,file,line);
  fflush(stderr);
  abort();
}

/*chech if the parameter has been found */
void check_parameter(const Flag_T flg,const char *const file, const int line)
{
  if (flg != FOUND)
  {
    fprintf(stderr,ERROR_MASSAGE
      "\"Parameter was not found.\"\n"
        "File: %s\nLine:%d\n",file,line);
    fflush(stderr);
    abort();
  }
  
}
