/*
// Alireza Rashti
// June 2018
*/

#include "text_tools.h"

/* strcmp case insensitive
// ->return value: 1 for success, 0 otherwise.
*/
int strcmp_i(const char *const s1, const char *const s2)
{
  if (!s2 || !s1) return 0;
  
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

/* strstr case insensitive
// ->return value: 1 for found, 0 otherwise.
*/
int strstr_i(const char *const s1, const char *const s2)
{
  if (!s2 || !s1) return 0;
  
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
    
  if (strstr(tmp1,tmp2))
    i = 1;
  else
    i = 0;
  
  free(tmp1);
  free(tmp2);
  
  return i;
}

/* duplicating a given string and INCLUDING '\0'
// ->return value: a pointer to string, or null if str is null
*/
char *dup_s(const char *const str)
{
  unsigned i,n;
  char *r = 0;
  
  if (!str) return 0;
  
  n = (unsigned)strlen(str)+1;
  
  r = malloc(n);
  pointerEr(r);
  
  for (i = 0; i < n-1; i++)
    r[i] = str[i];
  r[i] = '\0';
  
  return r;
}

/* 
// returning the first sub-string from str delimited by delimit
// and make savestr points to the the rest of str after delimit.
// NOTE: savestr must not be null.
// if str is 0 it uses savestr as the str.
// if str starts from delimit it'll be ignored and search begins from
// the next character.
//->return value: point to substring delimited by delimit, if nothing 
// is found point to NULL.
*/
char *tok_s(char *const str,const char delimit,char **const savestr)
{
  char *s = 0;
  char *ps = 0;
  unsigned l,i,f;
  Flag_T flg = NONE;

  assert(savestr);
  
  if (str != 0) s = str;
  else 		s = *savestr;
  
  /* return null if s is empty */
  if (!s) return 0;
  if (s[0] == '\0')  return 0;
  
  l = (unsigned)strlen(s);
  assert(s[l] == '\0');/* make sure this string end with null char */
  
  /* if start with delimit */
  if (s[0] == delimit) 
  {
    ++s;
    --l;
    flg = FOUND;
    *savestr = s;
    ps = s;
  }
  for (i = 0; i < l; i++)
  {
    if (s[i] == delimit)
    {
      /* make this sub-string empty of delimits */
      f = i;
      while (s[f] == delimit && s[f] != '\0' && f < l)
        s[f++] = '\0';
      
      ps = s;
      *savestr = &s[f];
      flg = FOUND;
      break;
    }
  }
  
  if (flg == NONE)
  {
    *savestr = 0;
    ps = s;
  }
    
  return ps;
}

/* taking out a portion of string delimited between d1 and d2.
// at their first occurrences.
// if str is 0 it uses savestr as the str.
// NOTE: savestr must not be null.
// ->return value: taken out substring and put the right leftover to save.
*/
char *sub_s(char *const str,const char d1,const char d2,char **const save)
{
  char *sub = 0;
  char *sub2 = 0;
  
  /* str = ...d1...d2...\0 */
  sub2 = tok_s(str,d2,save);/* sub2 = ...d1....(d2->\0)save */
  tok_s(sub2,d1,&sub); /* sub = (d1->\0+1)...(d2->\0) */
  
  return sub;
}
