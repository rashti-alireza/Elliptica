#include "core_lib.h"
#include "error_handling_lib.h"

/* global variables */
extern Parameter_T **parameters_global;

void update_parameter_double_format(const char *const lv, const double rv);
void update_parameter_integer(const char *const lv, const int rv);
void update_parameter_string(const char *const lv, const char *const rv);
void update_parameter_array(const char *const lv, const double *const rv,const unsigned n);
void read_input_file(const char *const path);
Parameter_T *get_parameter(const char *const par_name);
void add_parameter(const char *const lv, const char *const rv);
void add_parameter_string(const char *const lv, const char *const rv);
void add_parameter_double(const char *const lv, const double rv);
void add_parameter_array(const char *const lv, const double *const rv,const unsigned n);
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg);
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_double_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double *get_parameter_array_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
int make_parameters(const char *const path);
unsigned total_iterations_ip(void);
unsigned total_iterative_parameters_ip(void);
char *par_value_str_ip(const unsigned n);
char *par_name_ip(const unsigned n);
char *get_n_value_str_ip(const Parameter_T *const par,const unsigned n);
void update_iterative_parameter_ip(const unsigned iter);
unsigned total_iterative_parameters_ip(void);
unsigned total_iterations_ip(void);
static char *parse_multiplicity_of_iterative_parameter(const char *const rv);
void set_default_parameter(const char *const lhs,const char *const rhs);
void *alloc_parameter(Parameter_T ***const mem);
void free_parameter(const char *const par_name);
void free_given_parameter(Parameter_T *par);
void free_parameter_db(void);




