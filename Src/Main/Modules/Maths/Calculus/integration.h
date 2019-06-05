#include "core_lib.h"

void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
static double Composite_Simpson_1D(Integration_T *const I);



