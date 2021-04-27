/*
// Alireza Rashti
// June 2018
*/

#include "useful_functions.h"

/* helps you to find where little tests start*/
void test_start(const char *const file,const int line)
{
  printf("Test starts at\n"
    "File:%s\nLine:%d\n",file,line);
}

/* count the number of pointers which end to null 
// excluding the last one which is null
*/
Uint countf(void *const p)
{
  assert(p != 0);
  
  Uint c = 0;
  void **const pp = p;
  
  while (pp[c] != 0)
    c++;
    
  return c;
}

/* linear format to i component (row major order) */
Uint ijk_to_i_row_major_order(const Uint l, const Uint *const n)
{
  return l / (n[2]*n[1]);
}

/* linear format to j component (row major order) */
Uint ijk_to_j_row_major_order(const Uint l, const Uint *const n)
{
  Uint tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp / n[2];
}

/* linear format to k component (row major order) */
Uint ijk_to_k_row_major_order(const Uint l, const Uint *const n)
{
  Uint tmp;
  
  tmp = l % (n[2]*n[1]);
  return tmp % n[2];
}

/* changing text to enum for collocation 
// ->return value: collocation type, error if couldn't be found, 0 if null coll
*/
Collocation_T get_collocation(const char *const coll)
{
  Collocation_T c = UNDEFINED_COLLOCATION;
  
  if (!coll) return c;
  else if (strcmp_i(coll,"EquiSpaced")) c = EquiSpaced;
  else if (strcmp_i(coll,"Chebyshev_Extrema")) c = Chebyshev_Extrema;
  else if (strcmp_i(coll,"Chebyshev_Nodes")) c = Chebyshev_Nodes;
  else
    Errors("There is no such %s collocation.\n",coll);
    
  return c;
}

/* changing text to enum for basis 
// ->return value: basis type, error if couldn't be found, 0 if null basis
*/
Basis_T get_basis(const char *const basis)
{
  Basis_T b = UNDEFINED_BASIS;
  
  if (!basis) return b;
  else if (strcmp_i(basis,"Chebyshev_FirstKind")) 
    b = Chebyshev_Tn_BASIS;
  else if (strcmp_i(basis,"No_Basis")) 
    b = No_BASIS;
  else
    Errors("There is no such %s basis.\n",basis);
    
  return b;
}


/* find out if this point p located on an edge or not
// the algorithm is simple; if it happens at two or more interfaces,
// it means this point is on an edge and returns 1 otherwise 0
// note: it's only for points which are collocated
// ->return value: 1 if found 0 if not.
*/
int IsOnEdge(const Uint *const n,const Uint p)
{
  Uint i,j,k;
  int c;
  
  ijk_to_i_j_k(p,n,&i,&j,&k);
  
  c = 0;
  if (i == n[0]-1 || i == 0)  c++;
  if (j == n[1]-1 || j == 0)  c++;
  if (k == n[2]-1 || k == 0)  c++;
  
  if (c > 1)  return 1;
  
  return 0;
}

/* find out if this point p located on an face or not.
// the algorithm is simple; if it happens at one or more interfaces
// it means this point is on an face and returns number of face 
// otherwise 0.
// moreover, the found face f written like the example below:
// f[I_0] = 1, the point happens at face I_0,
// f[J_n1] = 0, the point won't happen at face J_n1 and etc.
// note: it's a general function and works for points which are not 
// collocated as well.
// ->return value: number of interface that found.
*/
int IsOnFace(const double *const x, const Patch_T *const patch,
             Uint *const f, const double precision_factor)
{
  int u,c;
  double X[3];
  
  for (u = 0; u < TOT_FACE; u++)
    f[u] = 0;
  
  c = 0;
  if(X_of_x_precision(X,x,patch,precision_factor))
  {
    for (u = 0; u < TOT_FACE; u++)
    {
      f[u] = check_interface(X,patch,u);
      
      if (f[u] == 1) c++;
    }
  }
  
  return c;
}

/* check if X is on specific interface.
// ->return value: 1 if yes, 0 otherwise.
*/
static Uint check_interface(const double *const X, const Patch_T *const patch, const int u)
{
  
  const Uint ind = patch->interface[u]->point[0]->ind;
  double res = EPS*root_square(3,X,0);
  double *Y;
  
  res = GRT(res,EPS) ? res: EPS;
  Y = patch->node[ind]->X;
  
  switch (u)
  {
    case I_0:
      if (LSSEQL(ABSd(Y[0]-X[0]),res)) return 1;
      break;
    case I_n0:
      if (LSSEQL(ABSd(Y[0]-X[0]),res)) return 1;
      break;
    case J_0:
      if (LSSEQL(ABSd(Y[1]-X[1]),res)) return 1;
      break;
    case J_n1:
      if (LSSEQL(ABSd(Y[1]-X[1]),res)) return 1;
      break;
    case K_0:
      if (LSSEQL(ABSd(Y[2]-X[2]),res)) return 1;
      break;
    case K_n2:
      if (LSSEQL(ABSd(Y[2]-X[2]),res)) return 1;
      break;
    default:
      Error0("No such interface was defined for this function.\n");
      
  }
  
  return 0;
}

/* given Cartesian value of a point, patch and f, it finds the closest
// node to this given point and return its index.
// ->return value: node index.
*/
Uint node_onFace(const double *const x, const Uint f,const Patch_T *const patch)
{
  const Interface_T *const face = patch->interface[f];
  Node_T **const node = patch->node;
  double nrm, s = DBL_MAX;
  Uint ind = UINT_MAX,i;
  
  FOR_ALL(i,face->point)
  {
    Uint id = face->point[i]->ind;
    nrm = root_square(3,x,node[id]->x);
    if (LSS(nrm,s))
    {
      s = nrm;
      ind = id;
    }
  }
  
  return ind; 
}

/* getting the subface which paired with sub
// ->return value: found paired subface.
*/
SubFace_T *get_paired_subface(const SubFace_T *const sub)
{
  const Patch_T *const patch = sub->patch->grid->patch[sub->adjPatch];
  const Interface_T *const face = patch->interface[sub->adjFace];
  
  assert(face);
  
  return face->subface[sub->adjsn];
  
}

/* ->return value: total number of nodes on the given grid */
Uint total_nodes_grid(const Grid_T *const grid)
{
  Uint p;
  Uint sum = 0;
  
  FOR_ALL_PATCHES(p,grid)
    sum += total_nodes_patch(grid->patch[p]);
  
  return sum;
}

/* ->return value: total number of nodes on the given patch */
Uint total_nodes_patch(const Patch_T *const patch)
{
  return patch->n[0]*patch->n[1]*patch->n[2];
}

/* find coord enum based in give str.
// ->return value: found Coord_T .
*/
Coord_T find_coord(const char *const coordsys)
{
  Coord_T coord = UNDEFINED_COORD;
  
  if(strcmp_i(coordsys,"Cartesian"))
    coord = Cartesian;
  else
    Errors("There is no such %s coordinates.\n",coordsys);  
    
  return coord;
}

/* generating a number number within final-initial.
// s determines weather srand be called or not.
// in fact, only when s is zero srand is called. and the reason is to avoid
// producing same random number in a loop.
// ->return value: random number in the double data type within the range.
*/
double random_double(const double initial,const double final,const Uint s)
{
  time_t t;
  
  /* Initializes random number generator */
  if (s == 0)
    srand((Uint) time(&t));
   
  return initial+((final-initial)/RAND_MAX)*rand();
}

/* hard copy subface s1 to subface s2 */
void copy_subface(SubFace_T *const s2,const SubFace_T *const s1)
{
  Uint i;

  s2->patch = s1->patch;
  s2->flags_str = dup_s(s1->flags_str);
  s2->sn = s1->sn;
  s2->adjsn = s1->adjsn;
  s2->np = s1->np;
  s2->id = calloc(s2->np,sizeof(*s2->id));
  IsNull(s2->id);
  s2->adjid = calloc(s2->np,sizeof(*s2->adjid));
  IsNull(s2->adjid);
  for (i = 0; i < s2->np; ++i)
    s2->id[i] = s1->id[i];
    
  if (s1->copy)
    for (i = 0; i < s2->np; ++i)
      s2->adjid[i] = s1->adjid[i];
      
  s2->face     = s1->face;
  s2->adjFace  = s1->adjFace;
  s2->adjPatch = s1->adjPatch;
  s2->df_dn  = s1->df_dn;
  s2->sameX  = s1->sameX;
  s2->sameY  = s1->sameY;
  s2->sameZ  = s1->sameZ;
  s2->touch  = s1->touch;
  s2->copy   = s1->copy;
  s2->exterF = s1->exterF;
  s2->outerB = s1->outerB;
  s2->innerB = s1->innerB;
}

/* map: points on each subface has an index for example id[i] = n.
// which says point with index i correspond to node n.
// this function gets n and return i mean invers(id[i]) = i;
// ->return value : invers(id[i])
*/
Uint subface_map_invers_id(const SubFace_T *const subface,const Uint n)
{
  Uint i;
  Uint s = UINT_MAX;
  
  for (i = 0; i < subface->np; ++i)
    if (subface->id[i] == n)
      return i;
  
  return s;  
}

/* given array s and its dimension, duplicate the array.
// ->return value: duplicated Uint array. */
Uint *dup_UINT(const Uint *const s,const Uint N)
{
  if (!s)
    Error0("The given array to be duplicated is Null!");
  
  Uint *dup = calloc(N,sizeof(*dup));
  IsNull(dup);
  Uint i;
  
  for (i = 0; i < N; ++i)
    dup[i] = s[i];
    
  return dup;
}

/* calculating the magnetude of d(X,Y,Z)/d(x,y,z) */
double max_Jacobian_dX_dx(Patch_T *const patch)
{
  double max = 0;
  double abs_j;
  const Uint nn = patch->nn;
  Uint l; 
  
  for (l = 0; l < nn; ++l)
  {
    abs_j = dq2_dq1(patch,_a_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_a_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_a_,_z_,l);
    if (abs_j > max)
      max = abs_j;
      
    abs_j = dq2_dq1(patch,_b_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_b_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_b_,_z_,l);
    if (abs_j > max)
      max = abs_j;  
      
    abs_j = dq2_dq1(patch,_c_,_x_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_c_,_y_,l);
    if (abs_j > max)
      max = abs_j;
    abs_j = dq2_dq1(patch,_c_,_z_,l);
    if (abs_j > max)
      max = abs_j;
      
  }
  
  return max;
}

/* calculating the error of spectral expansion using the last
// coefficient of the expansion for the given field f.
// ->return value: spectral expansion error */
double spectral_expansion_truncation_error(Field_T *const f)
{
  const double *const Cijk = make_coeffs_3d(f);
  const Uint *const n  = f->patch->n;
  
  return fabs(Cijk[i_j_k_to_ijk(n,n[0]-1,n[1]-1,n[2]-1)]);
}

/* go over all of the fields in the grid and print 
// the truncation error in this patch.*/
void print_spectral_expansion_truncation_error(Grid_T *const grid)
{
  pr_line_custom('=');
  printf("Truncation error at spectral expansion of the fields ...\n\n");

  const Uint np = grid->np;
  Uint p;
  
  UF_OpenMP(omp parallel for)
  for (p = 0; p < np; ++p)
  {
    Patch_T *patch = grid->patch[p];
    Uint nfld  = patch->nfld;
    double max_err = 0;
    const char *max_err_field_name = 0;
    const char *max_err_patch_name = 0;
    int len;
    Uint f;
    
    for (f = 0; f < nfld; ++f)
    {
      Field_T *field = patch->fields[f];
      
      /* skip residual fields */
      if (strstr(field->name,"_residual"))
        continue;
        
      double err  = spectral_expansion_truncation_error(field);
      
      if (err >= max_err)
      {
        max_err = err;
        max_err_field_name = field->name;
        max_err_patch_name = patch->name;
      }
    }
    if (!max_err_field_name || !max_err_patch_name)
      continue;
    len = 25-(int)strlen(max_err_field_name);
    printf("%s:\n",max_err_patch_name);
    printf("--> truncation error[%s] %*s %e\n",
              max_err_field_name,
              len,"= ",
              max_err);
  }/* end of for (p = 0; p < np; ++p) */
  
  printf("\nTruncation error at spectral expansion of the fields ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* calculating the maximum error of spectral derivative using the fact 
// that computer is using finite number of digits.
// f : given field
// o : given order of derivative
// ->return value: error in calculation = general idea is as follow:
// 1e-14*max(func)*max(Jacobian)^(order of derivative )*n*n^(2*order of derivative)*10  */
double spectral_derivative_max_error(const Field_T *const f,const Uint o)
{
  double e = 1e-15;
  double max_f,
         max_j;
  Uint max_n;
  const Uint *const n = f->patch->n;
  const char *der_par = PgetsEZ("Derivative_Method");
  
  if (strstr_i(der_par,"Spectral"))
  {
    max_f = L_inf(n[0]*n[1]*n[2],f->v);
    max_j = max_Jacobian_dX_dx(f->patch);
    max_n =  n[0] > n[1] ? n[0]  : n[1];
    max_n = max_n > n[2] ? max_n : n[2];
    e = 1e-14*max_n*pow(max_n,2*o);/* 1e-14 coming from: (machine precision=1e-15)*(all other unwarranted factors=10) */
    
    if (max_j > 1.0)/* if max_j is less than 1 it won't gonna make the error lesser! */
      e *= pow(max_j,o);
    if (max_f > 1.0)/* if max_f is less than 1 it won't gonna make the error lesser! */
      e *= max_f;
    if (e == 0.0) 
      e = 1e-15;
  } 
  else
    Error0(NO_JOB);
    
  return e;
}


/* ->return value: given stem name of a patch, 
// it finds the corresponding one in the given grid. */
Patch_T *GetPatch(const char *const stem,const Grid_T *const grid)
{
  Patch_T *retPatch = 0;
  char name[1000];
  Uint p;
  
  sprintf(name,PATCH_NAME_PRT_P_"%s",grid->gn,stem);
  
  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];
    
    if (strcmp_i(name,patch->name))
    {
      retPatch = patch;
      break;
    }
  }
  
  if(!retPatch)
    Errors("It could not find patch %s.\n",name);
    
  return retPatch;
}


/* print an array of double type with dimension n, for debuging purposes */
void dbprint(const double *v,const Uint n,const char *const desc)
{
  Uint i;
  
  pr_line();
  
  if (desc)/* description */
    printf("Debugging: %s\n",desc);
  else
    printf("Debugging:\n");
  
  for (i = 0; i < n; ++i)
    printf("data[%02u] = %+.15f\n",i,v[i]);
  
  pr_line();
}

/* issue a shell command using system */
void shell_command(const char *const cmd)
{
  int ret;
  printf("shell command:\n");
  printf("$  %s\n",cmd);
  fflush(stdout);
  ret = system(cmd);
  printf("=> returned (%d).\n",ret);
  fflush(stdout);
}

/* return value-> N*sizeof(double), using calloc*/
double *alloc_double(const Uint N)
{
  double *d;
  
  d = calloc(N,sizeof(*d));
  IsNull(d);
  
  return d;
}

/* return value-> M[R][C] double type memory using calloc */
double **alloc_2D_double(const long Uint R,const long Uint C)
{
  double **M;
  long Uint row;
  
  M = calloc(R,sizeof(*M));
  IsNull(M);
  
  for (row = 0; row < R; ++row)
  {
    M[row] = calloc(C,sizeof(*M[row]));
    IsNull(M[row]);
  }
  
  return M;
}

/* freeing 2 dimensions block of memory
// knowing the last column is pointing to null
*/
void free_2d(void *mem0)
{
  if (!mem0)
    return;
    
  int i;
  void **mem = mem0;
    
  for (i = 0; mem[i] != 0; i++)
    free(mem[i]);
    
  free(mem);
    
}

/* freeing 2 dimensions block of memory
// knowing the number of rows is c
*/
void free_2d_mem(void *mem0, const Uint long c)
{
  if (!mem0)
    return;
    
  Uint long i;
  void **mem = mem0;
  
  for (i = 0; i < c; i++)
    free(mem[i]);
    
  free(mem);
    
}

/* ->: memory usage, default Kb.
// The maximum resident set size used. That is, the maximum number of 
// physical memory that processes used simultaneously. */
double 
how_much_memory
  (
    const char *const unit/* gb,mb,kb */
  )
{
  struct rusage r_usage;
  double factor = 1;
  
  if (strcmp_i(unit,"gb"))
  {
    factor = 1E-6;
  }
  else if (strcmp_i(unit,"mb"))
  {
    factor = 1E-3;
  }
  else/* default */
  {
    factor = 1;
  }
    
  getrusage(RUSAGE_SELF,&r_usage);
  
  return ((double)r_usage.ru_maxrss*factor);
}

/* header and clock when the function called */
void header_and_clock(const char *const msg)
{
  const Uint Width = 74;/* number of rows in each line */
  const Uint lmsg  = (Uint)strlen(msg)+2/* "{ " */+5/* " ... " */;
  char clc[MAX_STR_LEN] = {'\0'};
  int d,h,m,s;
  
  /* clock */
  convert_clock(&d,&h,&m,&s);
  sprintf(clc,"[%02dd:%02dh:%02dm:%02ds]",d,h,m,s);
  
  if (Width>lmsg)/* right justification */
    printf("{ %s ... %*s\n",msg,Width-lmsg,clc);
  else
    printf("{ %s ... %s\n",msg,clc);
    
  fflush(stdout);
}

/* header and clock when the function ends */
void footer_and_clock(const char *const msg)
{
  const Uint Width = 74;/* number of rows in each line */
  const Uint lmsg  = (Uint)strlen(msg)+2/* "{ " */+5/* " ... " */;
  char clc[MAX_STR_LEN] = {'\0'};
  int d,h,m,s;
  
  /* clock */
  convert_clock(&d,&h,&m,&s);
  sprintf(clc,"[%02dd:%02dh:%02dm:%02ds]",d,h,m,s);
  
  if (Width>lmsg)/* right justification */
    printf("} %s :)) %*s\n",msg,Width-lmsg,clc);
  else
    printf("} %s :)) %s\n",msg,clc);
    
  fflush(stdout);
}

/* ->: field(X)|patch
// given field name, X and patch, interpolates the value of 
// the field at X on this patch */
double f_of_X(const char *const field_name,
              const double *const X/* patch coords */,
              Patch_T *const patch)
{
  double interp;
  Interpolation_T *interp_s = init_interpolation();
  Field_T *const F_field    = patch->fields[Ind(field_name)];
  
  interp_s->field = F_field;
  interp_s->X = X[0];
  interp_s->Y = X[1];
  interp_s->Z = X[2];
  interp_s->XYZ_dir_flag = 1;
  plan_interpolation(interp_s);
  interp = execute_interpolation(interp_s);
  free_interpolation(interp_s);
  
  return interp;
}

/* ->: max difference.
// calculate absolute difference between fields with stem1 and stem2. */
double diff_3x3_symmetric_fields(Grid_T *const grid,
                               const char *const stem1/* field1 */,
                               const char *const stem2/* field2 */,
                               const char *const rank/* [up/down] */,
                               const int pr_points/* print all points */)
{
  double max       = 0.;
  const int IsUp   = strcmp_i(rank,"up");
  const int IsDown = strcmp_i(rank,"down");
  Uint p;

  if (!IsUp && !IsDown)
    Errors("No such rank '%s' is defined.",rank);

  FOR_ALL_PATCHES(p,grid)
  {
    Patch_T *patch = grid->patch[p];

    if (IsUp)
    {
      READ_v_STEM(diff1_U2U2,stem1)
      READ_v_STEM(diff1_U1U2,stem1)
      READ_v_STEM(diff1_U1U1,stem1)
      READ_v_STEM(diff1_U0U2,stem1)
      READ_v_STEM(diff1_U0U1,stem1)
      READ_v_STEM(diff1_U0U0,stem1)
      
      READ_v_STEM(diff2_U2U2,stem2)
      READ_v_STEM(diff2_U1U2,stem2)
      READ_v_STEM(diff2_U1U1,stem2)
      READ_v_STEM(diff2_U0U2,stem2)
      READ_v_STEM(diff2_U0U1,stem2)
      READ_v_STEM(diff2_U0U0,stem2)
    
      FOR_ALL_ijk
      {
        CalcDiff(U2U2)
        CalcDiff(U1U2)
        CalcDiff(U1U1)
        CalcDiff(U0U2)
        CalcDiff(U0U1)
        CalcDiff(U0U0)
      }
    }
    else
    {
      READ_v_STEM(diff1_D2D2,stem1)
      READ_v_STEM(diff1_D1D2,stem1)
      READ_v_STEM(diff1_D1D1,stem1)
      READ_v_STEM(diff1_D0D2,stem1)
      READ_v_STEM(diff1_D0D1,stem1)
      READ_v_STEM(diff1_D0D0,stem1)
      
      READ_v_STEM(diff2_D2D2,stem2)
      READ_v_STEM(diff2_D1D2,stem2)
      READ_v_STEM(diff2_D1D1,stem2)
      READ_v_STEM(diff2_D0D2,stem2)
      READ_v_STEM(diff2_D0D1,stem2)
      READ_v_STEM(diff2_D0D0,stem2)
    
      FOR_ALL_ijk
      {
        CalcDiff(D2D2)
        CalcDiff(D1D2)
        CalcDiff(D1D1)
        CalcDiff(D0D2)
        CalcDiff(D0D1)
        CalcDiff(D0D0)
      }
    }
  }
  
  if (pr_points) 
  {
    printf(Pretty0"Linf{%s-%s} = %0.1e\n",stem1,stem2,max);
    fflush(stdout);
  }
  
  return max;
}

/* superimpose such that f = f1 + f2 + extra. */
void superimpose_simple(Grid_T *const grid,
                        const char *const f,
                        const char *const f1,
                        const char *const f2,
                        const double extra)
{
  OpenMP_Patch_Pragma(omp parallel for)
  FOR_ALL_p(grid->np)
  {
    Patch_T *patch = grid->patch[p];
    
    READ_v_STEM(F1,f1)
    READ_v_STEM(F2,f2)
    REALLOC_v_WRITE_v_STEM(F,f)
    
    FOR_ALL_ijk
      F[ijk] = F1[ijk]+F2[ijk]+extra;
  }
}


/* interpolate or copy the given comma separated fields 
// from old grid to new grid. 
// Note: when copy is used, it assumes the number of patches 
// and the resulotions are equal in both grids. */
void interpolate_fields_from_old_grid_to_new_grid
     (Grid_T *const ogrid/* old */,Grid_T *const ngrid/* new */,
     const char *const field_names/* comma separated field names */,
     const int copy/* if 1, only copy, if 0, only 3d interpolation */)
{
  FUNC_TIC
  
  /* some checks */
  if (!ogrid || !ngrid || !field_names)
  {
    FUNC_TOC
    printf(Pretty0"Empty argument(s)!\n");
    return;
  }
  
  char **fnames = read_separated_items_in_string(field_names,',');
  Uint f;
  
  /* if it is interpolation */
  if(copy == 0)
  {
    /* save points info for interpolation */
    struct Interp_S
    {
      double X[3];/* X of old patch */
      Patch_T *patch;/* old patch */
    }**pnts = 0;
    
    /* find all points and pertinent patches */
    pnts = calloc(ngrid->np,sizeof(*pnts));
    IsNull(pnts);
    OpenMP_Patch_Pragma(omp parallel for)
    FOR_ALL_p(ngrid->np)
    {
      Patch_T *patch = ngrid->patch[p];
      pnts[p] = calloc(patch->nn,sizeof(*pnts[p]));
      IsNull(pnts[p]);
      
      FOR_ALL_ijk
      {
        Patch_T *opatch = 
          x_in_which_patch(patch->node[ijk]->x,ogrid->patch,ogrid->np);

        /* approximate it */
        if (!opatch)
        {
          opatch = x_in_which_patch_force
            (patch->node[ijk]->x,ogrid->patch,ogrid->np,pnts[p][ijk].X);
          if(!opatch)
            Errors("I could not find a patch for %s.\nProbably a too "
                   "low resolution or strange configuration.",patch->name);
          pnts[p][ijk].patch = opatch;
          continue;
        }
        pnts[p][ijk].patch = opatch;
        X_of_x(pnts[p][ijk].X,patch->node[ijk]->x,opatch);
      }
    }
    
    /* make coeffs */
    f = 0;
    while(fnames[f])
    {
      OpenMP_Patch_Pragma(omp parallel for)
      FOR_ALL_p(ogrid->np)
      {
        Patch_T *patch = ogrid->patch[p];
        Field_T *fld   = patch->fields[Ind(fnames[f])];
        
        if (fld->v)
          make_coeffs_3d(fld);
        else
          Errors("There is no value for field '%s'.\n",fld->name);
      }
      f++;
    }
    
    /* interpolates from old to new */
    f = 0;
    while(fnames[f])
    {
      printf(Pretty0"Interpolating '%s' from grid%u to grid%u\n",
                                    fnames[f],ogrid->gn,ngrid->gn);
      fflush(stdout);
      
      OpenMP_Patch_Pragma(omp parallel for)
      FOR_ALL_p(ngrid->np)
      {
        Patch_T *patch = ngrid->patch[p];
        Field_T *fld   = patch->fields[Ind(fnames[f])];
        
        empty_field(fld);
        double *v = fld->v = alloc_double(patch->nn);
        
        FOR_ALL_ijk
        {
          v[ijk] = f_of_X(fnames[f],pnts[p][ijk].X,pnts[p][ijk].patch);
        }
      }
      f++;
    }
    
    /* free */
    if(pnts)
    FOR_ALL_p(ngrid->np)
      Free(pnts[p])
    Free(pnts);
  } 
  /* if copy */
  else if (copy == 1)
  {
    /* number of patches must be equal */
    if (ngrid->np != ogrid->np)
      Error0("Number of patches are different for copy!\n");

    f = 0;
    while(fnames[f])
    {
      printf(Pretty0"Copying '%s' from grid%u to grid%u\n",
             fnames[f],ogrid->gn,ngrid->gn);
      fflush(stdout);
      
      OpenMP_Patch_Pragma(omp parallel for)
      FOR_ALL_p(ngrid->np)
      {
        Patch_T *patch   = ngrid->patch[p];
        Field_T *fld_new = patch->fields[Ind(fnames[f])];
        Field_T *fld_old = ogrid->patch[p]->
                           fields[LookUpField_E(fnames[f],ogrid->patch[p])];
        assert(fld_old->v);
        empty_field(fld_new);
        double *v = fld_new->v = alloc_double(patch->nn);
        const char *name_p1 = strchr(patch->name,'_'),
                   *name_p2 = strchr(ogrid->patch[p]->name,'_');
        /* resolution must be equal */
        if(patch->nn != ogrid->patch[p]->nn)
          Error0("Resolution of patches are different for copy!\n");
        /* same patch must be used */  
        if (strcmp(name_p1,name_p2))
          Error0("Patch names are different for copy!\n");
        
        FOR_ALL_ijk
          v[ijk] = fld_old->v[ijk];
      }
      f++;
    }
  }
  else
    Error0(NO_OPTION);
  
  free_2d(fnames);
  
  FUNC_TOC
}
