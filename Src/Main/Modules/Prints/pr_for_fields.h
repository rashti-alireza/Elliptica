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
  


void pr_fields(Pr_Field_T *const pr);
Pr_Field_T *init_PrField(const Grid_T *const grid);
void free_PrField(Pr_Field_T *pr);
void pr_hdf5_silo(Pr_Field_T *const pr);
static void parse_parameter_3d(const char *const par,Pr_Field_T *const pr);
int print_fields_3D(const Grid_T *const grid,const int iteration,const char *const dir);
double print_fields_0D(const Grid_T *const grid,const int iteration,const char *const folder);
void print_fields_1D(const Grid_T *const grid,const int iteration,const char *const folder);
static Field_T **find_field_by_name_or_regex(const Patch_T *const patch,char **const fld_names,Uint *const Nfld);
static double map_to_patch_ref_interval(const double X,const Patch_T *const patch,const double *const min,const double *const max,const Uint dir, const int map_type);
static Uint find_closest_index(const double Xp,const Patch_T *const patch,const Uint dir);




