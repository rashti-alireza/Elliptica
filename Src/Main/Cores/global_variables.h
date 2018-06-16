#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilites.h"
#include <unistd.h>

#define MAX_ARR 400

/* global variables */
Parameter_T **parameters_global;// parameters
Project_T   **projects_global;// projects
char *path_global;// path of a directory where input file is
char *inputfile_name_global; //name of inputfile


static void make_path_global(void);
static void find_inputfile_name(char *const path);