/*
// Alireza Rashti
// August 2018
*/

#include "umfpack_method.h"

/* umfpack direct solver double interger kind 
// ->return value: EXIT_SUCCESS
*/
int direct_solver_umfpack_di(void *vp)
{
  Umfpack_T *const umf = vp;
  /* a is row by col matrix at a.x = b and in compressed column storage */
  const int row    = (int)umf->a->row;
  const int col    = (int)umf->a->col;
  int *const Ap    = umf->a->ccs->Ap;
  int *const Ai    = umf->a->ccs->Ai;
  double *const Ax = umf->a->ccs->Ax;
  double *const b  = umf->b;
  double *const x  = umf->x;
  void *Symbolic,*Numeric;
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  int status;
  double time1 = get_time_sec();
  
  if (!umf->a->ccs_f)
    Error0("The matrix a in a.x = b in umfpack must be in CCS format.\n");
  
  umfpack_di_defaults(Control);
  if (umf->refinement != DBL_MAX)/* refinment step */
    Control[UMFPACK_IRSTEP] = umf->refinement;
    
  status = umfpack_di_symbolic(row,col,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  umfpack_di_free_symbolic(&Symbolic);
  
  status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x,b,Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);

  /*freeing*/
  umfpack_di_free_numeric(&Numeric);
 
  umf->description[0] = '\0';
  sprintf(umf->description,
            "o.  Solver            = UMFPACK_di\n"
            "o.  Matrix Dimension  = %dx%d\n"
            "o.  Condition Number  = %e\n"
            "o.  Refinement_Step   = %d\n" 
            "o.  Solver Wall-Clock = %g(s)\n",
                                           row,col,
                                           1/Info[UMFPACK_RCOND],
                                           (int)Control[UMFPACK_IRSTEP],
                                           get_time_sec()-time1);
  
  return EXIT_SUCCESS;
}

/* umfpack direct solver double long interger kind 
// ->return value: EXIT_SUCCESS
*/
int direct_solver_umfpack_dl(void *vp)
{
  Umfpack_T *const umf = vp;
  /* a is row by col matrix at a.x = b and in compressed column storage */
  const int row    = (int)umf->a->row;
  const int col    = (int)umf->a->col;
  long *const Ap   = umf->a->ccs_long->Ap;
  long *const Ai   = umf->a->ccs_long->Ai;
  double *const Ax = umf->a->ccs_long->Ax;
  double *const b  = umf->b;
  double *const x  = umf->x;
  void *Symbolic,*Numeric;
  double Control[UMFPACK_CONTROL];
  double Info[UMFPACK_INFO];
  long status;
  double time1 = get_time_sec();
  const double Gb = 1.0E9;
  
  if (!umf->a->ccs_l_f)
    Error0("The matrix a in a.x = b in umfpack must be in CCS long format.\n");
  
  umfpack_dl_defaults(Control);
  if (umf->refinement != DBL_MAX)/* refinment step */
    Control[UMFPACK_IRSTEP] = umf->refinement;
  
  status = umfpack_dl_symbolic(row,col,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);
     
  status = umfpack_dl_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);
     
  umfpack_dl_free_symbolic(&Symbolic);

  status = umfpack_dl_solve(UMFPACK_A,Ap,Ai,Ax,x,b,Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);

  /*freeing*/
  umfpack_dl_free_numeric(&Numeric);

  umf->description[0] = '\0';
  sprintf(umf->description,
            "o.  Solver            = UMFPACK_dl\n"
            "o.  Matrix Dimension  = %dx%d\n"
            "o.  Condition Number  = %e\n"
            "o.  Sparsity          = %0.2f(%%)\n"
            "o.  Refinement Step   = %d\n"
            "o.  Peak Memory       = %g(Gbs)\n"
            "o.  Number of Flops   = %g\n"
            "o.  Number of Realloc = %g\n"
            "o.  Solver Wall-Clock = %g(s)\n",
             row,col,
             1/Info[UMFPACK_RCOND],
             (1-(double)Ap[col]/((double)row*col))*100.,
             (int)Control[UMFPACK_IRSTEP],
             Info[UMFPACK_PEAK_MEMORY]*Info[UMFPACK_SIZE_OF_UNIT]/Gb,
             Info[UMFPACK_FLOPS],
             Info[UMFPACK_NUMERIC_REALLOC],
             get_time_sec()-time1);

  return EXIT_SUCCESS;
}

/* umfpack was unsuccessful. it prints the reason and aborts */
void umfpack_error_di(const double *const Control,const int status,const char *const file,const int line)
{
  umfpack_di_report_status(Control,status);
  umfpack_failed(status,file,line);
}  

/* umfpack was unsuccessful. it prints the reason and aborts */
void umfpack_error_dl(const double *const Control,const long status,const char *const file,const int line)
{
  umfpack_dl_report_status(Control,status);
  umfpack_failed((int)status,file,line);
}  

/* umfpack direct solver double interger kind 
// used to solve ax = b for series of b's but same a.
// thus, it won't need to decompose matrix a each time.
// ->return value: EXIT_SUCCESS
*/
int direct_solver_series_umfpack_di(void *vp)
{
  Umfpack_T *const umf = vp;
  /* a is row by col matrix at a.x = b and in compressed column storage */
  const int row    = (int)umf->a->row;
  const int col    = (int)umf->a->col;
  int *const Ap    = umf->a->ccs->Ap;
  int *const Ai    = umf->a->ccs->Ai;
  double *const Ax = umf->a->ccs->Ax;
  double **const bs  = umf->bs;/* b's */
  double **const xs  = umf->xs;/* x's */
  const Uint ns = umf->ns;/* number of series to be solved */
  void *Symbolic,*Numeric;
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  int status;
  Uint i;
  double time1 = get_time_sec();
  
  if (!umf->a->ccs_f)
    Error0("The matrix a in a.x = b in umfpack must be in CCS format.\n");
  
  umfpack_di_defaults(Control);
  if (umf->refinement != DBL_MAX)/* refinment step */
    Control[UMFPACK_IRSTEP] = umf->refinement;
  
  status = umfpack_di_symbolic(row,col,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_di(Control,status,__FILE__,__LINE__);
     
  umfpack_di_free_symbolic(&Symbolic);

  /* solve for series ax[i]=b[i] */
  for (i = 0; i < ns; ++i)
  {
    status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,xs[i],bs[i],Numeric,Control,Info);
    if(status != UMFPACK_OK)
      umfpack_error_di(Control,status,__FILE__,__LINE__);
  }
  /*freeing*/
  umfpack_di_free_numeric(&Numeric);

  umf->description[0] = '\0';
  
  sprintf(umf->description,
            "o.  Series Solver     = UMFPACK_di\n"
            "o.  Matrix Dimension  = %dx%d\n"
            "o.  Condition Number  = %e\n"
            "o.  Refinement_Step   = %d\n" 
            "o.  Solver Wall-Clock = %g(s)\n",
                                           row,col,
                                           1/Info[UMFPACK_RCOND],
                                           (int)Control[UMFPACK_IRSTEP],
                                           get_time_sec()-time1);
                                           
  return EXIT_SUCCESS;
}

/* umfpack direct solver double long kind 
// used to solve ax = b for series of b's but same a.
// thus, it won't need to decompose matrix at each time.
// ->return value: EXIT_SUCCESS
*/
int direct_solver_series_umfpack_dl(void *vp)
{
  Umfpack_T *const umf = vp;
  /* a is row by col matrix at a.x = b and in compressed column storage */
  const int row     = (int)umf->a->row;
  const int col     = (int)umf->a->col;
  long *const Ap    = umf->a->ccs_long->Ap;
  long *const Ai    = umf->a->ccs_long->Ai;
  double *const Ax  = umf->a->ccs_long->Ax;
  double **const bs = umf->bs;/* b's */
  double **const xs = umf->xs;/* x's */
  const Uint ns = umf->ns;/* number of series to be solved */
  void *Symbolic,*Numeric;
  double Info[UMFPACK_INFO];
  double Control[UMFPACK_CONTROL];
  long status;
  Uint i;
  double time1 = get_time_sec();
  
  if (!umf->a->ccs_l_f)
    Error0("The matrix a in a.x = b in umfpack must be in CCS long format.\n");
  
  umfpack_dl_defaults(Control);
  if (umf->refinement != DBL_MAX)/* refinment step */
    Control[UMFPACK_IRSTEP] = umf->refinement;
  
  status = umfpack_dl_symbolic(row,col,Ap,Ai,Ax,&Symbolic,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);
     
  status = umfpack_dl_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);
  if(status != UMFPACK_OK)
    umfpack_error_dl(Control,status,__FILE__,__LINE__);
     
  umfpack_dl_free_symbolic(&Symbolic);

  /* solve for series ax[i]=b[i] */
  for (i = 0; i < ns; ++i)
  {
    status = umfpack_dl_solve(UMFPACK_A,Ap,Ai,Ax,xs[i],bs[i],Numeric,Control,Info);
    if(status != UMFPACK_OK)
      umfpack_error_dl(Control,status,__FILE__,__LINE__);
  }
  /*freeing*/
  umfpack_dl_free_numeric(&Numeric);
 
  umf->description[0] = '\0';
  
  sprintf(umf->description,
            "o.  Series Solver     = UMFPACK_dl\n"
            "o.  Matrix Dimension  = %dx%d\n"
            "o.  Condition Number  = %e\n"
            "o.  Refinement_Step   = %d\n" 
            "o.  Solver Wall-Clock = %g(s)\n",
                                           row,col,
                                           1/Info[UMFPACK_RCOND],
                                           (int)Control[UMFPACK_IRSTEP],
                                           get_time_sec()-time1);
                                           
  return EXIT_SUCCESS;
}

/* printing errors and reason of failure and then abort */
static void umfpack_failed(const int status,const char *const file,const int line)
{
  switch(status)
  {
    case UMFPACK_WARNING_singular_matrix:
      Error0("Matrix is singular. There are exact zeros on the diagonal of U.\n");
      break;
    case UMFPACK_WARNING_determinant_underflow:
      Error0("The determinant is nonzero, but smaller in\n"
          "magnitude than the smallest positive floating-point number.\n");
      break;
    case UMFPACK_WARNING_determinant_overflow:
      Error0("The determinant is larger in magnitude than\n"
        "the largest positive floating-point number (IEEE Inf).\n");
      break;
    case UMFPACK_ERROR_out_of_memory:
      Error0("Not enough memory. The ANSI C malloc or realloc\n"
        "routine failed.\n");
      break;
    case UMFPACK_ERROR_invalid_Numeric_object:
      Error0("Routines that take a Numeric object as input \n"
          "(or load it from a file) check this object and return this error code if it is invalid. This can\n"
          "be caused by a memory leak or overrun in your program, which can overwrite part of the\n"
          "Numeric object. It can also be caused by passing a Symbolic object by mistake, or some other\n"
          "pointer. If you try to factorize a matrix using one version of UMFPACK and then use the\n"
          "factors in another version, this error code will trigger as well. You cannot factor your matrix\n"
          "using version 4.0 and then solve with version 4.1, for example. 3 . You cannot use different\n"
          "precisions of the same version (real and complex, for example). It is possible for the Numeric\n"
          "object to be corrupted by your program in subtle ways that are not detectable by this quick\n"
          "check. In this case, you may see an UMFPACK_ERROR_different_pattern error code, or even\n"
          "an UMFPACK_ERROR_internal_error.\n");
      break;
    case UMFPACK_ERROR_invalid_Symbolic_object:
      Error0("Routines that take a Symbolic object as\n"
        "input (or load it from a file) check this object and return this error code if it is invalid.\n"
        "The causes of this error are analogous to the UMFPACK ERROR invalid Numeric object error\n"
        "described above.\n");
      break;
    case UMFPACK_ERROR_argument_missing:
      Error0("Some arguments of some are optional (you can pass a NULL pointer instead of an array).\n"
      "This error code occurs if you pass a NULL pointer when\n"
        "that argument is required to be present.\n");
      break;
    case UMFPACK_ERROR_n_nonpositive:
      Error0("The number of rows or columns of the matrix must be\n"
        "greater than zero.\n");
      break;
    case UMFPACK_ERROR_invalid_matrix:
      Error0("The matrix is invalid. For the column-oriented input,\n"
        "this error code will occur if the contents of Ap and/or Ai are invalid.\n"
        "Ap is an integer array of size n col+1. On input, it holds the “pointers” for the column form of\n"
        "the sparse matrix A. Column j of the matrix A is held in Ai [(Ap [j]) . . . (Ap [j+1]-1)].\n"
        "The first entry, Ap [0], must be zero, and Ap [j] ≤ Ap [j+1] must hold for all j in the\n"
        "range 0 to n col-1. The value nz = Ap [n col] is thus the total number of entries in the\n"
        "pattern of the matrix A. nz must be greater than or equal to zero.\n"
        "The nonzero pattern (row indices) for column j is stored in Ai [(Ap [j]) . . . (Ap [j+1]-1)].\n"
        "The row indices in a given column j must be in ascending order, and no duplicate row indices\n"
        "may be present. Row indices must be in the range 0 to n row-1 (the matrix is 0-based).\n"
        "Some routines take a triplet-form input, with arguments nz, Ti, and Tj. This error code is\n"
        "returned if nz is less than zero, if any row index in Ti is outside the range 0 to n col-1, or\n"
        "if any column index in Tj is outside the range 0 to n row-1.\n");
      break;
    case UMFPACK_ERROR_different_pattern:
      Error0("The most common cause of this error is that the\n"
        "pattern of the matrix has changed between the symbolic and numeric factorization. It can\n"
        "also occur if the Numeric or Symbolic object has been subtly corrupted by your program.\n");
      break;
    case UMFPACK_ERROR_invalid_system:
      Error0("The sys argument provided to one of the solve routines is invalid.\n");
      break;
    case UMFPACK_ERROR_invalid_permutation:
      Error0("The permutation vector provided as input is invalid.\n");
      break;
    case UMFPACK_ERROR_file_IO:
      Error0("This error code is returned by the routines that save and\n"
        "load the Numeric or Symbolic objects to/from a file, if a file I/O error has occurred. The file\n"
        "may not exist or may not be readable, you may be trying to create a file that you don’t have\n"
        "permission to create, or you may be out of disk space. The file you are trying to read might\n"
        "be the wrong one, and an earlier end-of-file condition would then result in this error.\n");
      break;
    case UMFPACK_ERROR_ordering_failed:
      Error0("The ordering method failed.\n");
      break;
    case UMFPACK_ERROR_internal_error:
      Error0("An internal error has occurred, of unknown cause.\n"
        "This is either a bug in UMFPACK, or the result of a memory overrun from your program. Try\n"
        "modifying the file AMD/Include/amd internal.h and adding the statement #undef NDEBUG,\n"
        "to enable the debugging mode. Recompile UMFPACK and rerun your program. A failed\n"
        "assertion might occur which can give you a better indication as to what is going wrong.\n"
        "Be aware that UMFPACK will be extraordinarily slow when running in debug mode. If all\n"
        "else fails, contact the developer (DrTimothyAldenDavis@gmail.com) with as many details as\n"
        "possible.\n");
      break;
    default:
      Error0("No such error defined for UMFPACK!\n");
  }
  UNUSED(file);
  UNUSED(line);
}

/* initializing the umfpack struct.
// ->return value: pointer to allocated umfpack struct */
Umfpack_T *init_umfpack(void)
{
  Umfpack_T *umf = calloc(1,sizeof(*umf));
  IsNull(umf);
  
  /* so this uses the default value of umfpack unless otherwise is set */
  umf->refinement = DBL_MAX;
  
  return umf;
}

/* free umfpack struct */
void free_umfpack(Umfpack_T *umf)
{
  if (!umf)
    return;
  free(umf);
}




