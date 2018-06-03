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
  char dir[MAX_ARR]     = {'\0'};
  char command[MAX_ARR] = {'\0'};
  char *path2 = 0;
  struct stat st = {0}; // status of files
  int i;
  
  i = 0;
  do 
  {
    sprintf(dir,"%s/%s_%d",path,name,i);
    i++;
  
  }while(stat(dir, &st) != -1);
  
  sprintf(command,"mkdir %s",dir);
  system(command);
  
  if (flg == YES)
  {
    path2 = strdup(dir);
  }
  
  return path2;
}