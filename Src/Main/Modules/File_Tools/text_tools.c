/*
// Alireza Rashti
// June 2018
*/

#include "text_tools.h"

/* strcmp case insensitive, it returns 1 for success otherwise 0 */
int strcmp_i(const char *const s1, const char *const s2)
{
  assert(s2 && s1);
  
  char *tmp1 = calloc(strlen(s1)+1,1);
  char *tmp2 = calloc(strlen(s2)+1,1);
  int i;
  
  i = 0;
  while(s1[i] != '\0')
  {
    tmp1[i] = (char)tolower(s1[i]);
    i++;
  }
  
  tmp1[i] = '\0';
  
  i = 0;
  while(s2[i] != '\0')
  {
    tmp2[i] = (char)tolower(s2[i]);
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

/* duplicating a given string and INCLUDING '\0'
// ->return value: a pointer to string
*/
char *dup_s(const char *const str)
{
  unsigned i,n;
  char *r = 0;
  
  assert(str);
  n = (unsigned)strlen(str)+1;
  
  r = malloc(n);
  pointerEr(r);
  
  for (i = 0; i < n-1; i++)
    r[i] = str[i];
  r[i] = '\0';
  
  return r;
}

/* 
// returning the first substring from str delimited by delimit
// and make savestr points to the the rest of str after delimit.
// if str is 0 it uses savestr as the str.
//->return value: point to substring delimited by delimit, if nothing 
// is found point to NULL.
*/
char *tok_s(char *const str,const char delimit,char **const savestr)
{
  char *s = 0;
  char *ps = 0;
  unsigned l,i;
  
  if (str != 0) s = str;
  else 		s = *savestr;
  
  l = (unsigned)strlen(s);
  
  for (i = 0; i < l; i++)
  {
    if (s[i] == delimit && s[i+1] != delimit && s[i+1] != '\0')
    {
      s[i] = '\0';
      ps = &s[i];
      *savestr = ps+1;
    }
  }
  
  return ps;
}
