/*
// Alireza Rashti
// June 2018
*/

#include "directory.h"

/* making a directory at the given path with the given name 
// and then returning the new folder path.
// ->return value: new folder path
*/
char *make_directory(const char *const path,const char *const name)
{
  char dir[MAX_ARR]     = {'\0'};
  char command[MAX_ARR] = {'\0'};
  char *folder = 0;
  struct stat st = {0};/* status of files */
  int i;
  
  i = 0;
  do 
  {
    sprintf(dir,"%s/%s_%d",path,name,i);
    i++;
  
  }while(stat(dir, &st) != -1);
  
  sprintf(command,"mkdir %s",dir);
  system(command);
  
  folder = dup_s(dir);
  
  return folder;
}