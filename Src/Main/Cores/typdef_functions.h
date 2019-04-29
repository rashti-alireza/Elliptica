/*
// Alireza Rashti
// June 2018
*/

/* *******************************************
// general purpose functions typedef
// *******************************************
*/

/* general function patch to void */
typedef void fFunc_PtoV_T (Patch_T *const patch);

/* general function grid to pointer to double */
typedef double *fFunc_Patch2Pdouble_T(Patch_T *const patch);

/* general function for variation of various kind of interpolation with 
// respect to the field. */
typedef double fdInterp_dfs_T(Patch_T *const patch,const double *const X,const unsigned df,const unsigned plane);

/* *******************************************
// general purpose structure typedef for functinos
// *******************************************
*/

/* patch to void struct */
typedef struct sFUNC_PtoV_T
{
  char *task;
  Coord_T coord;
  fFunc_PtoV_T *f;
}sFunc_PtoV_T;

/* grid to double struct */
typedef struct sFUNC_PATCH2PDOUBLE_T
{
  char *name;
  fFunc_Patch2Pdouble_T *func;
  unsigned flg: 1;/* used for different purposes */
}sFunc_Patch2Pdouble_T;

