#ifndef core_LIB_H
#define core_LIB_H


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
/* NOTE: elliptica_system_lib.h must include here */
#include "elliptica_system_lib.h"
#include "typdef_data.h"

/* some handy libraries */
#include "prints_lib.h"
#include "error_handling_lib.h"
#include "macros_lib.h"
#include "text_and_file_tools_lib.h"
#include "utilities_lib.h"

/* dealing with round off error */
#define ROUND_OFF_ERR (1E-12)
#define LSS(x,y) ((x) < (y)-ROUND_OFF_ERR)
#define GRT(x,y) ((x) > (y)+ROUND_OFF_ERR)
#define EQL(x,y) (((x) < (y)+ROUND_OFF_ERR) && ((x) > (y)-ROUND_OFF_ERR))
#define LSSEQL(x,y) (LSS(x,y) || EQL(x,y))
#define GRTEQL(x,y) (GRT(x,y) || EQL(x,y))

/* parameters */
void update_parameter_double_format(const char *const lv, const double rv,const int print_flg);
void update_parameter_integer(const char *const lv, const int rv);
void update_parameter_string(const char *const lv, const char *const rv);
void update_parameter_array(const char *const lv, const double *const rv,const Uint n);
void add_parameter(const char *const lv, const char *const rv);
void add_parameter_string(const char *const lv, const char *const rv);
void add_parameter_double(const char *const lv, const double rv,const int print_flg);
void add_parameter_array(const char *const lv, const double *const rv,const Uint n);
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg);
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_double_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double *get_parameter_array_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
Parameter_T *get_parameter(const char *const par_name);
Uint total_iterations_ip(void);
Uint total_iterative_parameters_ip(void);
char *par_value_str_ip(const Uint n);
char *par_name_ip(const Uint n);
char *get_n_value_str_ip(const Parameter_T *const par,const Uint n);
void update_iterative_parameter_ip(const Uint iter);
Uint total_iterative_parameters_ip(void);
Uint total_iterations_ip(void);
void set_default_parameter(const char *const lhs,const char *const rhs);
void *alloc_parameter(Parameter_T ***const mem);
void free_parameter(const char *const par_name);
void free_given_parameter(Parameter_T *par);
void free_parameter_db(void);

int 
update_iteration_params
  (
   const Uint main_loop_iter,
   const char *const prefix/* parameter prefix */,
   const char *const dir_name_format/* eg: "BHNS_%s_%ux%ux%u" */
  );

char *par_name_rule(const char *const s);

#endif




