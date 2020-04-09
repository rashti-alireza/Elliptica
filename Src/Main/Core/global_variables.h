#include "core_lib.h"
#include "error_handling_lib.h"
#include "utilities_lib.h"
#include <unistd.h>

#define MAX_ARR 400
#define EXTENSION ".in"

/* global variables */
Grid_T **grids_global;
Parameter_T **parameters_global;/* parameters */
Project_T   **projects_global;/* projects */
time_t initial_time_global;/* initial time abc starts */

int init_global_variables(const char *const path);
static void find_relative_root_path(const char *const path);
static void find_inputfile_name(const char *const path);
