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
#include "print_lib.h"
#include "macros_lib.h"
#include "text_tools_lib.h"

#define ROUND_OFF_ERR 1E-10
#define LSS(x,y) (x < y-ROUND_OFF_ERR)
#define GRT(x,y) (x > y+ROUND_OFF_ERR)
#define EQL(x,y) (x < y+ROUND_OFF_ERR && x > y-ROUND_OFF_ERR)
#define LSSEQL(x,y) (LSS(x,y) || EQL(x,y))
#define GRTEQL(x,y) (GRT(x,y) || EQL(x,y))

/* parameters */
int    get_parameter_value_I(const char *const par_name,Flag_T *const flg);
double get_parameter_value_D(const char *const par_name,Flag_T *const flg);
char  *get_parameter_value_S(const char *const par_name,Flag_T *const flg);
/* fields */
Field_T *init_field_3d(const char *const name,Grid_T *const grid);
void add_field(Field_T *const f,Grid_T *const grid);
Field_T *get_field_S(const char *const name,Grid_T *const grid);
double *make_coeffs(Field_T *const f);
double *make_coeffs_inverse(Field_T *const f);
