#include "core_lib.h"/* it includes prints_lib.h too */
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include "macros_lib.h"
#include "maths_general_lib.h"
#include "manifold_lib.h"
#include "fields_lib.h"
#include "pr_for_fields_structs.h"
#include  <unistd.h>

#define STR_LEN (999)

/* openmp */
#define PR_1D_OpenMP(x) _Pragma ( #x )
#define PR_2D_OpenMP(x) _Pragma ( #x )

/* compactify repeated lines, dir = (dir)ection */
#define FWRITE_1D_VALUES(dir) \
  fprintf(file,"%0.15f %0.15f %0.15f %0.15f",\
               patch->node[ijk]->X[dir],\
               patch->node[ijk]->x[0],\
               patch->node[ijk]->x[1],\
               patch->node[ijk]->x[2]);\
  for (f = 0; f < Nfld; ++f){\
    fprintf(file," %0.15f",fields[f]->v[ijk]);}\
  fprintf(file,"\n");

/* compactify repeated lines, dir = (dir)ection */
#define FWRITE_2D_VALUES(dir1,dir2) \
  fprintf(file,"%0.15f %0.15f %0.15f %0.15f %0.15f",\
               patch->node[ijk]->X[dir1],\
               patch->node[ijk]->X[dir2],\
               patch->node[ijk]->x[0],\
               patch->node[ijk]->x[1],\
               patch->node[ijk]->x[2]);\
  for (f = 0; f < Nfld; ++f){\
    fprintf(file," %0.15f",fields[f]->v[ijk]);}\
  fprintf(file,"\n");

/* general 1D header, dir = (dir)ection */
#define FWRITE_1D_HEADER(dir) \
  if (dir == 0)\
  {\
   fprintf(file,"line_%s: t in [%0.2f,%0.2f] -->"\
                " t*I + (Y[%u]=%0.2f)*J + (Z[%u]=%0.2f)*K\n",\
                  line->strv,patch->min[dir],patch->max[dir],\
                  J,Y,K,Z);\
  }\
  else if (dir == 1)\
  {\
   fprintf(file,"line_%s: t in [%0.2f,%0.2f] -->"\
                " (X[%u]=%0.2f)*I + t*J + (Z[%u]=%0.2f)*K\n",\
                  line->strv,patch->min[dir],patch->max[dir],\
                  I,X,K,Z);\
  }\
  else if (dir == 2)\
  {\
   fprintf(file,"line_%s: t in [%0.2f,%0.2f] -->"\
                " (X[%u]=%0.2f)*I + (Y[%u]=%0.2f)*J + t*K\n",\
                  line->strv,patch->min[dir],patch->max[dir],\
                  I,X,J,Y);\
  }\
  else Error0(NO_OPTION);
  

/* general 2D header, pln = {XYplane = 0, XZplane = 1, YZplane = 2} */
#define FWRITE_2D_HEADER(pln) \
  if (pln == 0)\
  {\
   fprintf(file,"plane_%s: t1 x t2 in [%0.2f,%0.2f]x[%0.2f,%0.2f] -->"\
                " t1*I + t2*J + (Z[%u]=%0.2f)*K\n",\
                  plane->strv,patch->min[0],patch->max[0],\
                  patch->min[1],patch->max[1],K,Z);\
  }\
  else if (pln == 1)\
  {\
   fprintf(file,"plane_%s: t1 x t2 in [%0.2f,%0.2f]x[%0.2f,%0.2f] -->"\
                " t1*I + (Y[%u]=%0.2f)*J + t2*K\n",\
                  plane->strv,patch->min[0],patch->max[0],\
                  patch->min[2],patch->max[2],J,Y);\
  }\
  else if (pln == 2)\
  {\
     fprintf(file,"plane_%s: t1 x t2 in [%0.2f,%0.2f]x[%0.2f,%0.2f] -->"\
                " (X[%u]=%0.2f)*I + t1*J + t2*K\n",\
                  plane->strv,patch->min[1],patch->max[1],\
                  patch->min[2],patch->max[2],I,X);\
  }\
  else Error0(NO_OPTION);
  


void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
void pr_hdf5_silo(Pr_Field_T *const pr);
static void parse_parameter_3d(const char *const par,Pr_Field_T *const pr);
int print_fields_3D(const Grid_T *const grid,const int iteration,const char *const dir);
double print_fields_0D(const Grid_T *const grid,const int iteration,const char *const folder);
void print_fields_1D(const Grid_T *const grid,const int iteration,const char *const folder);
void print_fields_2D(const Grid_T *const grid,const int iteration,const char *const folder);
static Field_T **find_field_by_name_or_regex(const Patch_T *const patch,char **const fld_names,Uint *const Nfld);
static double map_to_patch_ref_interval(const double X,const Patch_T *const patch,const double *const min,const double *const max,const Uint dir, const int map_type);
static Uint find_closest_index(const double Xp,const Patch_T *const patch,const Uint dir);




