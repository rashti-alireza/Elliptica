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

void *get_parameter_value(char *const par_name,Flag_T kind, double *value);
