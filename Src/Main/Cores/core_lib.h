#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "typdef_data.h"
#include "typdef_functions.h"
#include "prints_lib.h"
#include "error_handling_lib.h"
#include "macros_lib.h"
#include "text_tools_lib.h"

#define ROUND_OFF_ERR 1E-12
#define LSS(x,y) (x < y-ROUND_OFF_ERR)
#define GRT(x,y) (x > y+ROUND_OFF_ERR)
#define EQL(x,y) (x < y+ROUND_OFF_ERR && x > y-ROUND_OFF_ERR)
#define LSSEQL(x,y) (LSS(x,y) || EQL(x,y))
#define GRTEQL(x,y) (GRT(x,y) || EQL(x,y))

/* parameters */
void add_parameter(const char *const lv, const char *const rv);
void add_parameter_string(const char *const lv, const char *const rv);
void add_parameter_double(const char *const lv, const double rv);
void add_parameter_array(const char *const lv, const double *const rv,const unsigned n);
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg);
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_double_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double *get_parameter_array_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
Parameter_T *get_parameter(const char *const par_name);
unsigned total_iterations_ip(void);
unsigned total_iterative_parameters_ip(void);
char *par_value_str_ip(const unsigned n);
char *par_name_ip(const unsigned n);
char *get_n_value_str_ip(const Parameter_T *const par,const unsigned n);
void update_iterative_parameter_ip(const unsigned iter);
unsigned total_iterative_parameters_ip(void);
unsigned total_iterations_ip(void);

/* fields */
Field_T *add_field(const char *const name,const char *attribute,Patch_T *const patch,const Flag_T alloc_flg);
void remove_field(Field_T *f);
void add_attribute(Field_T *const fld,const char *const attribute);
int LookUpField(const char *const name,const Patch_T *const patch);
double *make_coeffs_1d(Field_T *const f,const unsigned dir);
double *make_coeffs_2d(Field_T *const f,const unsigned dir1,const unsigned dir2);
double *make_coeffs_3d(Field_T *const f);
void enable_fields(Grid_T *const grid);

/* grid */
Patch_T make_temp_patch(const Patch_T *const patch);
void free_temp_patch(Patch_T *const patch);

