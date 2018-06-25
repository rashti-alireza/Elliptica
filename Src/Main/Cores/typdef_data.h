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
  INUSE,
  FOUND
}Flag_T;

/* collocation */
typedef enum COLLOCATION_T
{
  EquiSpaced,
  Chebyshev_Zero
}Collocation_T;

/* coordinate system */
typedef enum COORD_T
{
  Cartesian,
  CubedSphere
}Coord_T;

/* print flags */
typedef enum PRINT_T
{
  PRINT_PARAMETER,
  PRINT_COORDS
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

/* node*/
typedef struct NODE_T
{
  double cart[3];// for Cartesian value x,y,z
  double *curv;// for general curvilinear value a,b,c
}Node_T;

/* point */
typedef struct POINT_T
{
  int ind;// linear index
  double N[3];// normal vector
  struct PATCH_T *patch;// refers to its patch
  int adjPatch;// adjacent patch used in interpolation
  struct POINT_T *adjPoint;// adjacent point
  unsigned int face   : 5;// the interface in which this point located
  unsigned int touch  : 1;// touch state 1, overlap state 0
  unsigned int copy   : 1;// copy state 1, interpolation state 0
  unsigned int exterF : 1;// external interface 1, internal 0
                          // extenal means it can reach other interface
  unsigned int outerB : 1;// if it is outer boundary of grid 
                          // and needs boundary condition
  unsigned int innerB : 1;// if it is inner boundary of grid 
                          // and needs boundary condition                 
}Point_T;

/* face */
typedef struct Interface_T
{
  int np;// number of points in this structure
  Point_T **point;// points on the interface
}Interface_T;

/* patch */
typedef struct PATCH_T
{
  struct GIRD_T *grid;// referes to its grid
  char *name;// box name
  char *coordsys;// coord sys used in this patch
  Collocation_T collocation;// type of collocation in this patch
  int n[3];// number of point in each direction
  double c[3];// center
  double s[3];// size like length, width and height
  double min[3];// minimum of each direction like x_min = min[0]
  double max[3];// maximum of each direction like b_max = max[1]
  Node_T **node;// node info
  Interface_T *interface;// interface info 
}Patch_T;

/* grid */
typedef struct GIRD_T
{
  char *kind;// type of grid which refers how we cover the grid
  Flag_T status;// INUSE or READY
  Patch_T **patch;// covering patch
  
  //Field_T 
}Grid_T;
