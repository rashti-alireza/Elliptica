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
  YES,
  NO,
  NONE,
  READY,
  INUSED,
  FOUND
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

/* functions */
typedef int ProjFunc(void);

/* projects */
typedef struct PROJECT_T
{
  char *name;// name of project
  char *des;// description
  ProjFunc *func;// project function
}Project_T;

/* *******************************************
// grid and related: 
// *******************************************
*/

/* nodes*/
typedef struct NODE_T
{
  double cart[3];// for Cartesian value x,y,z
  double *curv;// for general curvilinear value a,b,c
}Node_T;

/* patch */
typedef struct PATCH_T
{
  char *name;// box name
  char *coordsys;// coord sys used in this patch
  char *collocation;// type of collocation in this patch
  int n[3];// number of point in each direction
  double c[3];// center
  double s[3];// size like length, width and height
  double min[3];// minimum of each direction like x_min = min[0]
  double max[3];// maximum of each direction like b_max = max[1]
  Node_T **node;// nodes info
}Patch_T;

/* grid */
typedef struct GIRD_T
{
  char *kind;// type of grid wich refers how we cover the grid
  Flag_T status;// INUSE or READY
  Patch_T **patch;// covering patch
  
  //Field_T 
}Grid_T;
