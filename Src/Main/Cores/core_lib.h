#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "typdef_data.h"
#include "typdef_functions.h"
#include "print_lib.h"
#include "macros_lib.h"

int    get_parameter_value_I(char *const par_name,Flag_T *flg);
double get_parameter_value_D(char *const par_name,Flag_T *flg);
char  *get_parameter_value_S(char *const par_name,Flag_T *flg);
