/* struct fot integration */
typedef struct INTEGRATION_T
{
 const char *type;
 double err;/* expected error this method */
 struct
 {
  double a,b;/* integral limits */
  double *f;/* integrand */
  unsigned n;/* odd positive integer  */
 }Composite_Simpson_1D[1];
 double (*integration_func)(struct INTEGRATION_T *const I);/* function that integrates */
}Integration_T;

double *Partial_Derivative(Field_T *const f,const char *task);
void free_integration(Integration_T *I);
void plan_integration(Integration_T *const I);
double execute_integration(Integration_T *const I);
Integration_T *init_integration(void);
int integration_tests(Grid_T *const grid);

