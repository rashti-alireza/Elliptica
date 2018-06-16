/*
// Alireza Rashti
// May 2018
*/

/* *******************************************
// enum:
// *******************************************
*/

/* common flags*/
typedef enum FLAG_T
{
  NUMERIC,
  LITERAL,
  YES,
  NONE
}Flag_T;

/* print flags */
typedef enum PRINT_T
{
  PRINT_PARAMETER
}Print_T;

/* *******************************************
// parameter:
// *******************************************
*/

/* parameters */
typedef struct PARAMETER_T
{
  /* syntax is expected to be lv = rv */
  char *lv;//letf value
  char *rv;// right value
}Parameter_T;

/* *******************************************
// project:
// *******************************************
*/

/* functions*/
typedef int ProjFunc(void);

/* projects */
typedef struct PROJECT_T
{
  char *name;// name of project
  char *des;// description
  ProjFunc *func;// project function
}Project_T;

/* *******************************************
// grid:
// *******************************************
*/