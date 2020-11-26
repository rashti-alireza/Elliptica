#include "sys_header.h"


typedef void fAdjustment_t (Physics_T *const phys);


int sys_tune_ADM_momenta(Physics_T *const phys);
static fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust);
static void parse_adjust_parameter(const char *const par,char *adjust[3]);
static void Px_ADM_is0_by_y_CM(Physics_T *const phys);
static void Py_ADM_is0_by_x_CM(Physics_T *const phys);




