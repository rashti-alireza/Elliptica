#ifndef fields_LIB_H
#define fields_LIB_H
#include "elliptica_system_lib.h"


/* forward declaration structures */
struct GRID_T;
struct PATCH_T;

/* Field structure */
typedef struct FIELD_T
{
  char *name;/* name of field */
  double *v;/* values on each node on patch */
  double *v2;/* other values. if this field has two kinds of value:
             // e.g. fields in spectral expansion needs both 
             // coefficients of expansion and values of field on each node.
             // so for field with expansion this v2 refers to coefficents
             // value. */
  double *v3;/* other values. for example in Ylm expansion we need 
             // both real and imaginary coeffs, one goes with v2,
             // the other goes with v3; */
  char *attr;/* attributes of fields like its dimension 
             // or other essential info.
             */
  char *info;/* each field might need more info or attributes 
             // which will save here and they are temporary.
             // e.g. the info about the coeffs of field. this info 
             // is dynamic.
             */
  struct PATCH_T *patch;/* refers to its patch which this field defined */
  void *v_ptr;/* refers to various quantities based on need */
}Field_T;


Field_T *add_field(const char *const name,const char *attribute,struct PATCH_T *const patch,const Flag_T alloc_flg);
void remove_field(Field_T *f);
void add_attribute(Field_T *const fld,const char *const attribute);
int LookUpField(const char *const name,const struct PATCH_T *const patch);
int LookUpField_E(const char *const name,const struct PATCH_T *const patch);
double *make_coeffs_1d(Field_T *const f,const Uint dir);
double *make_coeffs_2d(Field_T *const f,const Uint dir1,const Uint dir2);
double *make_coeffs_3d(Field_T *const f);
void enable_fields(struct GRID_T *const grid);
void free_v2(Field_T *f);
void free_v3(Field_T *f);
void free_attr(Field_T *f);
void free_v(Field_T *f);
void free_info(Field_T *f);
void free_field(Field_T *fld);
void empty_field(Field_T *fld);
int free_coeffs(Field_T *fld);
void remove_field_regex(Patch_T *const patch,const char *const regex);
Uint *find_field_index_regex(const Patch_T *const patch,const char *const regex,Uint *const Nm);
void remove_aux_fields(struct GRID_T *const grid,const char *const aux_names);
void add_aux_fields(struct GRID_T *const grid,const char *const aux_names);

#endif






