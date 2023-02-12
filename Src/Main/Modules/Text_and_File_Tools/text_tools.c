/*
// Alireza Rashti
// June 2018
*/

#include "text_tools.h"

#define MAX_STR_LEN 400

/* strcmp case insensitive
// note: if either of them is null 0 is returned.
// ->return value: 1 for success, 0 otherwise.
*/
int strcmp_i(const char *const s1/* haystack */, 
             const char *const s2/* needle */)
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
// note: if either of them is null 0 is returned.
// ->return value: 1 for found, 0 otherwise.
*/
int strstr_i(const char *const s1/* haystack */, 
             const char *const s2/* needle */)
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
  Uint i,n;
  char *r = 0;
  
  if (!str) return 0;
  
  n = (Uint)strlen(str)+1;
  
  r = calloc(n,sizeof(*r));
  IsNull(r);
  
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
  Uint l,i,f;
  Flag_T flg = NONE;

  assert(savestr);
  
  if (str != 0) s = str;
  else 		s = *savestr;
  
  /* return null if s is empty */
  if (!s) return 0;
  if (s[0] == '\0')  return 0;
  
  l = (Uint)strlen(s);
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

/* getting a format and string, find out if the string is in compliance
// with the format. in case there are numbers of sub-strings with the
// supposedly same format, it ONLY checks the first sub-string.
// format is like: "{?(?)|?}" so it checks to see the string has all of the
// delimits specified in the format. NOTE: ONLY DELIMITS MUST BE SPECIFIED.
// WORDS are TO BE WRITTEN BY ?.
// ->return value: 1 if in compliance, 0 otherwise.
*/
int check_format_s(const char *str,const char *const format)
{
  int r;
  char (*delimits)[2];
  char trim[MAX_STR_LEN]  = {'\0'}/* triming the format */,
       trim2[MAX_STR_LEN] = {'\0'}/* triming the str */;
  char *subs1,*subs2;/* make sure you check the first substring */
  Uint i,n;
  
  if(!str)
    return 0;
  delimits = 0;
  
  n = 0;
  i = 0;
  while (format[i] != '\0')
  {
    if (format[i] != '?')
    {
      delimits = realloc(delimits,(n+1)*sizeof(*delimits));
      IsNull(delimits);
      delimits[n][0] = format[i];
      delimits[n][1] = '\0';
      trim[n] = format[i];
      trim[n+1] = '\0';
      n++;
    }
    i++;
  }
  
  subs1 = dup_s(str);
  subs2 = strchr(subs1,delimits[n-1][0]);
  if (!subs2)
    return 0;
    
  subs2[1] = '\0';
  str = subs1;
  
  r = 1;
  for (i = 0; i < n; i++)
  {
    str = strchr(str,delimits[i][0]);
    if (str == 0)
    {
      r = 0;
      break;
    }
    else
    {
      trim2[i] = str[0];
      trim2[i+1] = '\0';
      str++;
    }
  }
  
  if (strcmp(trim,trim2))
  {
    r = 0;
  }
  
  if (delimits)
    free(delimits);
  free(subs1);  
  
  return r;
}

/* given a lists (heystack) and number (N) of names,
// it returns the index of the matched needle.(case sensitive)
// ->return value: i if heystach[i] = needle, otherwise UINT_MAX.
*/
Uint find_index_string(char **const heystack,const Uint N,const char *const needle)
{
  Uint i,j = UINT_MAX;
  
  for (i = 0; i < N; ++i)
    if (!strcmp(heystack[i],needle))
    {
      j = i;
      break;
    }
  
  return j;
}

/* it checks if the pattern matches in the str 
// by using regular expression libraries.
// it's case sensitive, and new line is NOT treated as separator of new string.
// note: regex_pattern must be in posix format, i.e. [:digit:], [:alpha:] etc.
// ex: regex_pattern = ".+(left|right)_NS_around.+"
//     => matches: _left_NS_around_
// ->return value: 1 if matches, 0 otherwise */
int regex_search(const char *const regex_pattern,const char *const str)
{
  regex_t regex;
  int ret_regex;
  int ret = 0;
  
  ret_regex = regcomp(&regex,regex_pattern,REG_EXTENDED);
  if (ret_regex)
    Error0("Regular expression Faild. It could be wrong pattern or memory failure.\n");
  
  ret_regex = regexec(&regex,str,0,0,0);
  
  if (!ret_regex)
    ret = 1;

  regfree(&regex);
      
  return ret;
}

// it finds and returns a copy of the first match.
// if there is no match, then return NULL.
// it's case sensitive, and new line is NOT treated as separator of new string.
// note: regex_pattern must be in posix format, i.e. [:digit:], [:alpha:] etc.
// ->return value: return the first match, NULL otherwise */
char *regex_find(const char *const regex_pattern,const char *const str)
{
  regex_t regex;
  char *ret = 0;
  const Uint n_matches = 1;/* number of matches */
  regmatch_t match[n_matches];
  int status;
  int len = 0;
  int i;
  
  status = regcomp(&regex,regex_pattern,REG_EXTENDED);
  if (status)
    Error0("Regular expression Faild. It could be wrong pattern or memory failure.\n");
  
  status = regexec(&regex,str,n_matches,match,0);
  if (status)/* if no match is found */
  {
    regfree(&regex);
    return 0;
  }
  
  regfree(&regex);
  
  len = match[0].rm_eo/* The offset in string of the end of the substring.  */
      - match[0].rm_so/* The offset in string of the beginning of a substring. */
      ;
  ret = calloc((Uint long)len+1,1);
  IsNull(ret);
  
  for (i = 0; i < len; ++i)
    ret[i] = str[match[0].rm_so+i];
  ret[len] = '\0';
  
  return ret;
}

/* parsing a string contains items separated with delimiter and return
// all items separately. note: the last pointer is null so one can count
// the number of items. note: it conserves the order from the left of give string.
// ->return value: array of pointers to items,the last pointer is null. */
char **read_separated_items_in_string(const char *const string,const char delimiter)
{
  /* if null return null */
  if (!string)
    return 0;
    
  char *str = dup_s(string);/* str = f1,f2,... */
  char **items = 0;
  Uint ni  = 0;/* number of items */
  char *tok,*save = 0;

  tok = tok_s(str,delimiter,&save);/* tok = f1 */
  while(tok)
  {
    items = realloc(items,(ni+2)*sizeof(*items));
    IsNull(items);
    items[ni]   = dup_s(tok);
    items[ni+1] = 0;
    tok = tok_s(0,delimiter,&save);/* tok = f2 */
    ni++;
  }
  free(str);
  
  return items;
}

/* ->: 1 if found match and replace, 0 otherwise.
// replace the matches match of given string orig by repl
// and write into save.
// NOTE: save must have enough memory.  */
int regex_replace(const char *const orig/* original */,
                  const char *const regex_pattern/* regex pattern */,
                  const char *const repl/* replace by this piece */,
                  char *const save/* write the result in save  */)
{
  regex_t regex;
  char *original = dup_s(orig);
  const char *c  = original;
  const Uint n_matches = 1;/* number of matches */
  regmatch_t match[n_matches];
  int status;
    
  status = regcomp(&regex,regex_pattern,REG_EXTENDED);
  if (status)
    Error0("Regular expression Faild. It could be wrong pattern or memory failure.\n");

  status = regexec(&regex,c,n_matches,match,0);
  
  if (status)/* if no match is found write original to save. */
  {
    regfree(&regex);
    sprintf(save,"%s",original);
    Free(original);
    return 0;
  }
  else/*  there is at least one match */
  {
    Uint i = 0;
    int len = match[0].rm_eo/* The offset in string of the end of the substring.  */
             -match[0].rm_so;/* The offset in string of the beginning of a substring. */
    while(!status)
    {
      const char *o  = c;
      
      /* write the leading */
      while (c != &o[match[0].rm_so])
      {
        save[i] = *c;
        c++;
        i++;
      }
      Uint j = 0;
      /* replace */
      while (repl[j] != '\0')
        save[i++] = repl[j++];
      
      c += len;/* skip match */
      status = regexec(&regex,c,n_matches,match,0);
      
      len = match[0].rm_eo-match[0].rm_so;
      /* if regex is "^" */
      if (len == 0) 
        status = 1;/* => breaks the loop */
    }
    /* write the leading */
    while (*c != '\0')
    {
      save[i] = *c;
      i++;
      c++;
    }
    save[i] = '\0';
    regfree(&regex);
  }
  
  Free(original);
  return 1;
}


