double *Partial_Derivative(Field_T *const f,const char *task);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
int integration_tests(Grid_T *const grid);

