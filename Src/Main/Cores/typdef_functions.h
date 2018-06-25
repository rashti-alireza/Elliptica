/*
// Alireza Rashti
// June 2018
*/

/* *******************************************
// general purpose functions and related struct
// *******************************************
*/

/* general function patch to void */
typedef void fFunc_PtoV_T(Patch_T *patch);

/* patch to void struct */
typedef struct sFUNC_PtoV_T
{
  char *task;
  Coord_T coord;
  fFunc_PtoV_T *f;
}sFunc_PtoV_T;

