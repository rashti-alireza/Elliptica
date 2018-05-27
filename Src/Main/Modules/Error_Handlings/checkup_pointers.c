/*
// Alireza Rashti
// May 2018
*/

#include "checkup_pointers.h"

/* checking up if the allocation of the pointer is failed;
// if it is failed then exit.
*/
void checkup_pointer(void *const p, char *file, int line)
{
  if (p == NULL)
  {
    fprintf(stderr,"!!!!!!ERROR!!!!!!:\n"
      "\"Pointer allocation was failed.\n\""
        "File: %s\nLine:%d\n");
    printf("!!!!!!EXIT due to ERROR!!!!!!\n");
    abort();
  }
}

/* bad input */
void bad_input_error(char *file, int line)
{
    fprintf(stderr,"!!!!!!ERROR!!!!!!:\n"
      "\"There is no such input for the function.\n\""
        "File: %s\nLine:%d\n");
    printf("!!!!!!EXIT due to ERROR!!!!!!\n");
    abort();
}