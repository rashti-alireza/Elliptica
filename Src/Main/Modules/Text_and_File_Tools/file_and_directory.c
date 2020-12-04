/*
// Alireza Rashti
// June 2018
*/

#include "file_and_directory.h"

/* making a directory at the given path with the given name 
// and then returning the new folder path.
// ->return value: new folder path
*/
char *make_directory(const char *const path,const char *const name)
{
  char dir[MAX_ARR]     = {'\0'};
  char command[MAX_ARR2] = {'\0'};
  char *folder = 0;
  struct stat st = {0};/* status of files */
  int i;
  
  i = 0;
  do 
  {
    sprintf(dir,"%s/%s"FOLDER_AFFIX,path,name,i);
    i++;
  }while(!stat(dir, &st));
  
  sprintf(command,"mkdir %s",dir);
  shell_command(command);
  
  folder = dup_s(dir);
  
  return folder;
}

/* making folder in "top_directory". 
// note: allocate memory for path.
// ->return value: made folder directory
*/
char *make_folder(const char *const folder)
{
  const char *path_par;
  char *path;
  path_par = Pgets("top_directory");
  path = make_directory(path_par,folder);
  
  return path;
}

/* open the given folder name if already exists otherwise make one.
// ->return value: path of specified folder.
*/
char *open_folder(const char *const folder)
{
  struct stat st = {0};/* status of files */
  const char *path_par = Pgets("top_directory");
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

/* open file and check if it is successful.
// ->return: opened file */
void *fopen_and_check(const char *const file_path,const char *const mode,const char *const file_dbg, const int line_dbg)
{
  FILE *file = fopen(file_path,mode);
  
  if (!file)
    abort_error_string("Could not open file '%s'\n",file_path,file_dbg,line_dbg);
  
  return file;
}

