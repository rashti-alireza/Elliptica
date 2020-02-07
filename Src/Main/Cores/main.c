/*
// Alireza Rashti
// May 2018
*/

#include "main.h"

int main(int argn, char **argv)
{
  pr_logo();
  /* if no parameter file is given */
  if (argn == 1)
  {
    pr_info();
    return EXIT_SUCCESS;
  }
  
  /* run the project */
  const char *proj_name;
  Project_T *proj;
    
  init_global_variables(argv[argn-1]);/* initiating global variables */
  make_parameters(argv[argn-1]);/* reading input file 
                                // and making parameters */
  
  create_db_projects();/* add all of the desired projects 
                       // to the project data base */
  
  proj_name = PgetsEZ("Project");
  proj      = get_project(proj_name);
  
  /* check if the parameter file is correct and this project exists */
  if (!proj)
    abortEr_s("There is no such %s project!\n",proj_name);
  
  execute_project(proj);
  
  /* free the data base of project */
  free_db_project();
  
  return EXIT_SUCCESS;
}

/* print logo */
static void pr_logo(void)
{
  const char *const logo =
  
" Welcome to\n\n\n"  
" EEEEEEEEEEEE  ll  ll\n"
" EE            ll  ll\n"             
" EE            ll  ll  O                       O\n"              
" EE            ll  ll                   tt\n"
" EEEEEEEEEEEE  ll  ll  ii  pppppppp. tttttttt  ii    .ccccc   .aaaaa.\n"
" EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     a.\n"
" EE            ll  ll  ii  pp       p.  tt     ii  .c       .a       a.\n"
" EE            ll  ll  ii  pp      p.   tt     ii   .c       .a     .a\n"
" EEEEEEEEEEEE  ll  ll  ii  pppppppp.    tt     ii    .ccccc   .aaaaa.aa\n"
"                           pp                                          aa\n"  
"                           pp\n"
"                           pp\n"
"                           pp\n"
;

  pr_line_custom('*');
  printf("%s\n",logo);
  printf("\nBug reports and feedback to: rashti.alireza@gmail.com\n");
  pr_line_custom('*');
}

/* some info about the Elliptica */
static void pr_info(void)
{
  printf("\n"
    "  Elliptica is an infrastructure to solve boundary value problems in\n"
    "general relativity to construct initial data for compact objects.\n\n"
    "usage:\n$ ./elliptica parameter_file.par\n\n");

 fflush(stdout);
}
