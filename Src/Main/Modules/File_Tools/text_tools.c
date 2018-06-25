/*
// Alireza Rashti
// June 2018
*/

#include "text_tools.h"

/* strcmp case insensitive, it returns 1 for success otherwise 0 */
int strcmp_i(char *s1, char *s2)
{
  assert(s2 && s1);
  
  char *tmp1 = calloc(strlen(s1)+1,1);
  char *tmp2 = calloc(strlen(s2)+1,1);
  int i;
  
  i = 0;
  while(s1[i] != '\0')
  {
    tmp1[i] = tolower(s1[i]);
    i++;
  }
  
  tmp1[i] = '\0';
  
  i = 0;
  while(s2[i] != '\0')
  {
    tmp2[i] = tolower(s2[i]);
    i++;
  }
  tmp2[i] = '\0';
    
  if (!strcmp(tmp1,tmp2))
    i = 1;
  else
    i = 0;
  
  free(tmp1);
  free(tmp2);
  
  return i;
}
