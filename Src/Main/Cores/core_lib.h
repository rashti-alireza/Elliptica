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
#include "text_tools_lib.h"

#define ROUND_OFF_ERR 1E-10
#define LSS(x,y) (x < y-ROUND_OFF_ERR)
#define GRT(x,y) (x > y+ROUND_OFF_ERR)
#define EQL(x,y) (x < y+ROUND_OFF_ERR && x > y-ROUND_OFF_ERR)
#define LSSEQL(x,y) (LSS(x,y) || EQL(x,y))
#define GRTEQL(x,y) (GRT(x,y) || EQL(x,y))

int    get_parameter_value_I(char *const par_name,Flag_T *flg);
double get_parameter_value_D(char *const par_name,Flag_T *flg);
char  *get_parameter_value_S(char *const par_name,Flag_T *flg);
