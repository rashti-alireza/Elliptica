#include "core_lib.h"
#include "error_handling_lib.h"

#define PAR_LEN  (300)
#define STR_SIZE1 (99)
#define STR_SIZE2 (999)

/* print param related */
#define PAR_WIDTH_PR  (29)
#define PAR_FORMAT_PR "%-*s = %+5.4e,|%+5.4e|,(%+08.3f%%)"

/* string format for param's name */
#define PAR_NAME_RULE(Xname) par_rule1_uppercase_lowercase(Xname)

/* prefix the given string str with the given prefix pre and 
// return string prepar.
// NOTE: prepar must be defined. */
#define PrefixIt(pre,str) (sprintf(prepar,"%s%s",pre,str) ? prepar : 0)


/* global variables */
extern Parameter_T **parameters_global;

void update_parameter_double_format(const char *const lv, const double rv,const int print_flg);
void update_parameter_integer(const char *const lv, const int rv);
void update_parameter_string(const char *const lv, const char *const rv);
void update_parameter_array(const char *const lv, const double *const rv,const Uint n);
void read_input_file(const char *const path);
Parameter_T *get_parameter(const char *const par_name);
void add_parameter(const char *const lv, const char *const rv);
void add_parameter_string(const char *const lv, const char *const rv);
void add_parameter_double(const char *const lv, const double rv,const int print_flg);
void add_parameter_array(const char *const lv, const double *const rv,const Uint n);
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg);
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double get_parameter_double_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
double *get_parameter_array_format(const char *const par_name,const char *const file, const int line,const Flag_T flg);
int make_parameters(const char *const path);
Uint total_iterations_ip(void);
Uint total_iterative_parameters_ip(void);
char *par_value_str_ip(const Uint n);
char *par_name_ip(const Uint n);
char *get_n_value_str_ip(const Parameter_T *const par,const Uint n);
void update_iterative_parameter_ip(const Uint iter);
Uint total_iterative_parameters_ip(void);
Uint total_iterations_ip(void);
static char *parse_multiplicity_of_iterative_parameter(const char *const rv);
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

static char *par_rule1_uppercase_lowercase(const char *const s) __attribute__ ((unused));
static char *par_rule2_lowercase(const char *const s) __attribute__ ((unused));
char *par_name_rule(const char *const s);




