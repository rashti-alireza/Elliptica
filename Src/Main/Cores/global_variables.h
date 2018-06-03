#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilites.h"
#include <unistd.h>

#define MAX_ARR 400

/* global variables */
Parameter_T **global_parameter;// parameters
char *global_path;// path of a directory where input file is
char *global_inputfile_name; //name of inputfile


static void make_global_path(void);
static void find_inputfile_name(char *const path);