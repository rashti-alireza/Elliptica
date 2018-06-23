#include "core_lib.h"
#include "error_handling_lib.h"
#include "memory_managing_lib.h"

/* global variables */
extern Parameter_T **parameters_global;
extern char *path_global;

void read_input_file(char *const path);
Parameter_T *get_parameter(char *const par_name);
void add_parameter(char *lv, char *rv);
char *make_directory(char *path,char *name,Flag_T flg);
void set_default_parameter(void);

