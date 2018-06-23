#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include <unistd.h>

#define MAX_ARR 400
#define EXTENSION ".in"

/* global variables */
Grid_T **grids_global;
Parameter_T **parameters_global;// parameters
Project_T   **projects_global;// projects
time_t initial_time_global;// initial time abc starts
char *path_global;// path of a directory where input file is
char *inputfile_name_global; //name of inputfile


static void make_path_global(void);
static void find_inputfile_name(char *const path);