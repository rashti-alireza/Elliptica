#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"

/* global variables */
extern Parameter_T **parameters_global;
extern char *path_global;

void read_input_file(const char *const path);
Parameter_T *get_parameter(const char *const par_name);
void add_parameter(const char *const lv, const char *const rv);
void make_directory(char **const path,const char *const name,const Flag_T flg);
void set_default_parameter(void);
int get_parameter_value_I(const char *const par_name,Flag_T *const flg);
double get_parameter_value_D(const char *const par_name,Flag_T *const flg);
char *get_parameter_value_S(const char *const par_name,Flag_T *const flg);
int make_parameters(const char *const path);

