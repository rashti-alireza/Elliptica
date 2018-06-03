/*
// Alireza Rashti
// June 2018
*/

#include "directory.h"

/* making a directory with given path and 
// returning the path of made directory if flg = YES
*/
char *make_directory(char *path,char *name,Flag_T flg)
{
  char dir[1000]={'\0'};
  char *path2 = 0;
  
  sprintf(dir,"mkdir %s\\%s",path,name);
  system(dir);
  
  if (flg == YES)
  {
    path2 = strdup(dir);
  }
  
  return path2;
}