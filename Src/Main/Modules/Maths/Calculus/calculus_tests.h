#include "core_lib.h"
#include "memory_managing_lib.h"
#include "maths_general_lib.h"
#include "maths_calculus_lib.h"
#include "maths_analytic_lib.h"

#define DO 1
#define NOT_DO 0

int integration_tests(Grid_T *const grid);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
static int csr_1d(Grid_T *const grid);
static int GQ_ChebExtrema(Grid_T *const grid);
static int GQ_Lobatto(Grid_T *const grid);
