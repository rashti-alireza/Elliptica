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
    sprintf(dir,"%s/%s"FOLDER_AFFIX,path,name,i);
    i++;
  }while(!stat(dir, &st));
  
  printf("shell command:\n$ mkdir %s\n\n",dir);
  fflush(stdout);
  sprintf(command,"mkdir %s",dir);
  system(command);
  
  folder = dup_s(dir);
  
  return folder;
}

/* making folder in "output_directory_path". 
// note: allocate memory for path.
// ->return value: made folder directory
*/
char *make_folder(const char *const folder)
{
  const char *path_par;
  char *path;
  path_par = Pgets("output_directory_path");
  path = make_directory(path_par,folder);
  
  return path;
}

/* open the given folder name if already exists otherwise make one.
// ->return value: path of specified folder.
*/
char *open_folder(const char *const folder)
{
  struct stat st = {0};/* status of files */
  const char *path_par = Pgets("output_directory_path");
  char dir[MAX_ARR];
  char *ret = 0;
  int i;
  
  i = 0;
  do 
  {
    sprintf(dir,"%s/%s_%04d",path_par,folder,i);
    i++;
  
  }while(stat(dir, &st) != -1);
  
  /* if the folder already exists */
  if (i > 1)
  {
    sprintf(dir,"%s/%s_%04d",path_par,folder,i-2);
    ret = dup_s(dir);
  }
  /* if not exists */
  else
  {
    ret = make_folder(folder);
  }
  
  return ret;
}
