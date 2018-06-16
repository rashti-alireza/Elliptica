/*
// Alireza Rashti
// May 2018
*/

#include "checkup_pointers.h"

/* checking up if the allocation of the pointer is failed;
// if it is failed then exit.
*/
void checkup_pointer_error(void *const p, char *file, int line)
{
  if (p == NULL)
  {
    fprintf(stderr,ERROR_MASSAGE
      "\"Pointer allocation was failed.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf(ERROR_MASSAGE_ABORT);
    abort();
  }
}

/* bad input */
void bad_input_error(char *file, int line)
{
    fprintf(stderr,ERROR_MASSAGE
      "\"There is no such input for the function.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf(ERROR_MASSAGE_ABORT);
    abort();
}

/* null path directory */
void null_path_error(void *path,char *file, int line)
{
  if (path == NULL)
  {
    fprintf(stderr,ERROR_MASSAGE
      "\"The directory path is null.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf(ERROR_MASSAGE_ABORT);
    abort();
  }
}

/* general purpose error */
void abort_error(char *massage,char *file, int line)
{
    fprintf(stderr,ERROR_MASSAGE"%s"
        "File: %s\nLine:%d\n",massage,file,line);
    printf(ERROR_MASSAGE_ABORT);
    abort();
}

/* general purpose error with more detail string version */
void abort_error_string(char *massage,char *detail,char *file, int line)
{
  char msg[400] = {'\0'};
  
  sprintf(msg,massage,detail);
  
  fprintf(stderr,ERROR_MASSAGE"%s"
      "File: %s\nLine:%d\n",msg,file,line);
  printf(ERROR_MASSAGE_ABORT);
  abort();
}
