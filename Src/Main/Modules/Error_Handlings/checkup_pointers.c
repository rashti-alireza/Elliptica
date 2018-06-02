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
    fprintf(stderr,"\t!!!!!!ERROR!!!!!!\n"
      "\"Pointer allocation was failed.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf("\t!!!!!!EXIT due to ERROR!!!!!!\n");
    abort();
  }
}

/* bad input */
void bad_input_error(char *file, int line)
{
    fprintf(stderr,"\t!!!!!!ERROR!!!!!!\n"
      "\"There is no such input for the function.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf("\t!!!!!!EXIT due to ERROR!!!!!!\n");
    abort();
}

/* null path directory */
void null_path_error(char *path,char *file, int line)
{
  if (path == NULL)
  {
    fprintf(stderr,"\t!!!!!!ERROR!!!!!!\n"
      "\"The directory path is null.\"\n"
        "File: %s\nLine:%d\n",file,line);
    printf("\t!!!!!!EXIT due to ERROR!!!!!!\n");
    abort();
  }
}