#ifndef maths_matrix_solvers
#define maths_matrix_solvers

struct MATRIX_T;

/* umfpack direct solver */
typedef struct UMFPACK_T
{
  char description[9000];
  struct MATRIX_T *a;/* a in a.x = b */
  double *b;/* in ax=b */
  double *x;/* in ax=b */
  double **xs;/* x for series solving Ax[i]=b[i], i runs 0,...,ns */
  double **bs;/* b for series solving Ax[i]=b[i], i runs 0,...,ns */
  unsigned ns;/* number of series */
}UmfPack_T;

int direct_solver_umfpack_di(void *vp);
int direct_solver_umfpack_dl(void *vp);
int direct_solver_series_umfpack_di(void *vp);
int direct_solver_series_umfpack_dl(void *vp);
void solver_tests(void);

#endif

