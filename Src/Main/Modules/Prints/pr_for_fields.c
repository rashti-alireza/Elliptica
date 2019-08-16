/*
// Alireza Rashti
// July 2018
*/

#include "pr_for_fields.h"

/* DeLimits */
#define DL_OC '{'
#define DL_CC '}'
#define DL_OP '('
#define DL_CP ')'
#define DL_C ','

/* Print Fields:
//
// usage examples:
// ===============
// # parameter that is determined in input file is like follows:
// print_fields_4d = yes,Format:HDF5,{(V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2)}
// Pr_Field_T *pr  = init_PrField(grid);
// pr->folder      = "folder_path";
// pr->par         = "print_fields_4d";
// pr_fields(pr);
// free_PrField(pr);
*/

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and print the result in the specified folder over ALL GRID.
*/
void pr_fields(Pr_Field_T *const pr)
{
  const char *par4 = GetParameterS(pr->par);/* print_fields_4d */
  const char *par3 = GetParameterS(pr->par);/* print_fields_3d */
  const char *par2 = GetParameterS(pr->par);/* print_fields_2d */
  
  /* 4d prints, i.e. field versus a,b,c */
  if (strstr_i(par4,"yes"))
  {
    read_parameter_4d(par4,pr);
    if (strstr_i(par4,"Format:HDF5"))
    {
      pr_fields_on_grid_HDF5_4d(pr);
    }
    else
      abortEr("The print format for 4d print is "
        "not recognized in the input-file.\n");
    
  }/* end of if (strstr_i(par4,"yes")) */
  /* 3d prints, i.e. ?*/
  else if (strstr_i(par3,"yes"))
  {
    abortEr(INCOMPLETE_FUNC);
  }/* end of if (strstr_i(par3,"yes")) */
  /* 2d prints, i.e. ?*/
  else if (strstr_i(par2,"yes"))
  {
    abortEr(INCOMPLETE_FUNC);
  }/* end of if (strstr_i(par2,"yes")) */
  
}

/* freeing info struct */
static void free_info_s(Pr_Field_T *const pr)
{
  struct Info_S *info = pr->group;
  unsigned i;
  
  for (i = 0; i < pr->ng; ++i)
  {
    if (info[i].vec_flg)
    {
      free(info[i].comp[0]);
      free(info[i].comp[1]);
      free(info[i].comp[2]);
    }
    else
      free(info[i].field);
    
  }
  free(info);

}

/* printing fields with HDF5 format using silo library */
static void pr_fields_on_grid_HDF5_4d(Pr_Field_T *const pr)
{
  struct Info_S *const pr_info = pr->group;
  const unsigned npr = pr->ng;
  unsigned pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    DBfile *dbfile_Cart = 0;
    DBfile *dbfile_Curv = 0;
    unsigned i_pr;
   
    /* printing Cartesian mesh */
    dbfile_Cart = make_structured_mesh_3d_Cartesian(pr,patch);
    /* printing Curvilinear mesh */
    dbfile_Curv = make_structured_mesh_3d_Curvilinear(pr,patch);
    
    /* printing field on on the meshes */
    for (i_pr = 0; i_pr < npr; ++i_pr)
    {
      pr->file = dbfile_Cart;
      pr->file2 = dbfile_Curv;
      pr->vptr = &pr_info[i_pr];
      /* if it is vector to be printed */
      if (pr_info[i_pr].vec_flg)
        pr_vector_on_structured_mesh_3d_silo(pr);
      
      /* printing scalar field on mesh */
      else
        pr_scalar_on_structured_mesh_3d_silo(pr);
    }
    
    if (DBClose(dbfile_Cart) == -1)
      abortEr("Silo library failed to close the file.\n");
    
    if (DBClose(dbfile_Curv) == -1)
      abortEr("Silo library failed to close the file.\n");
  }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
}

/* printing 3D mesh with Cartesian value 
// ->return value: the file containing the mesh.
*/
static void *make_structured_mesh_3d_Cartesian(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  float *x = 0,*y = 0,*z = 0;
  char *label[3];
  prepare_node_structured_mesh_3d_silo("Cartesian",patch,&x,&y,&z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s.Cart%04d.silo",
    pr->folder,patch->name,pr->cycle);
  dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
    "3D mesh in Cartesian values with HDF5 format using silo library",
    DB_HDF5);
  pointerEr(dbfile);

  /* setting up some options */
  label[0] = dup_s("x");
  label[1] = dup_s("y");
  label[2] = dup_s("z");
  DBoptlist *opt = DBMakeOptlist(5);
  DBAddOption(opt,DBOPT_XLABEL,(void *)label[0]);
  DBAddOption(opt,DBOPT_YLABEL,(void *)label[1]);
  DBAddOption(opt,DBOPT_ZLABEL,(void *)label[2]);
  DBAddOption(opt,DBOPT_DTIME,&pr->time);
  DBAddOption(opt,DBOPT_CYCLE,&pr->cycle);
  
  /* setting up pr */
  pr->a = x;
  pr->b = y;
  pr->c = z;
  pr->patch = patch;
  pr->opt_patch = opt;
  pr->file = dbfile;
  /* printing mesh */
  pr_structured_mesh_3d_silo(pr);
  free_nodes_silo(x,y,z);
  DBFreeOptlist(opt);
  free(label[0]);
  free(label[1]);
  free(label[2]);
  
  return dbfile;
}

/* printing 3D Curvilinear mesh 
// ->return value: the file containing the mesh.
*/
static void *make_structured_mesh_3d_Curvilinear(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  float *x = 0,*y = 0,*z = 0;
  char *label[3];
  prepare_node_structured_mesh_3d_silo("Curvilinear",patch,&x,&y,&z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s.Curv%04d.silo",
    pr->folder,patch->name,pr->cycle);
  dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
    "3D Curvilinear mesh with HDF5 format using silo library",
    DB_HDF5);
  pointerEr(dbfile);

  /* setting up some options */
  label[0] = dup_s("a");
  label[1] = dup_s("b");
  label[2] = dup_s("c");
  DBoptlist *opt = DBMakeOptlist(5);
  DBAddOption(opt,DBOPT_XLABEL,(void *)label[0]);
  DBAddOption(opt,DBOPT_YLABEL,(void *)label[1]);
  DBAddOption(opt,DBOPT_ZLABEL,(void *)label[2]);
  DBAddOption(opt,DBOPT_DTIME,&pr->time);
  DBAddOption(opt,DBOPT_CYCLE,&pr->cycle);
  
  /* setting up pr */
  pr->a = x;
  pr->b = y;
  pr->c = z;
  pr->patch = patch;
  pr->opt_patch = opt;
  pr->file = dbfile;
  /* printing mesh */
  pr_structured_mesh_3d_silo(pr);
  free_nodes_silo(x,y,z);
  DBFreeOptlist(opt);
  free(label[0]);
  free(label[1]);
  free(label[2]);

  return dbfile;
}


/* printing scalar on structured 3d mesh with node centered format.
// note data for silo library must be written in column format,
// i.e. 3d to 1d map is i+n0*(j+n1*k).
*/
static void pr_scalar_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  assert(pr->file);
  
  DBfile *const dbfile = pr->file; 
  const Patch_T *const patch = pr->patch;
  struct Info_S *subg = pr->vptr;
  float *data;
  double *v;
  const int nnodes = (int)total_nodes_patch(patch);
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  const int v_ind = Ind(subg->field);
  int DB_ret;
  unsigned i,j,k;
  int l;
  
  if (v_ind < 0)
    abortEr_s("There is no such field \"%s\" among the fields!\n",subg->field);
  
  data = malloc((long unsigned)nnodes*sizeof(*data));
  pointerEr(data);
  
  v = patch->pool[v_ind]->v;
  for (l = 0 ; l < nnodes; ++l)
  {
    IJK((unsigned)l,patch->n,&i,&j,&k);
    data[i+(unsigned)dims[0]*(j+(unsigned)dims[1]*k)] = (float)v[l];
  }
   
  DB_ret = DBPutQuadvar1(dbfile,subg->field,pr->patch->name,
    data,dims,ndims,0,0,DB_FLOAT,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
  /* if there is another file that the field needs to be printed */
  if (pr->file2)
    DB_ret = DBPutQuadvar1(pr->file2,subg->field,pr->patch->name,
      data,dims,ndims,0,0,DB_FLOAT,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
  free(data);
}

/* printing vector field on structured 3d mesh with node centered format.
// note data for silo library must be written in cloumn format,
// i.e. 3d to 1d map is i+n0*(j+n1*k).
*/
static void pr_vector_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  assert(pr->file);
  
  DBfile *const dbfile = pr->file; 
  const Patch_T *const patch = pr->patch;
  struct Info_S *subg = pr->vptr;
  float *comp[3];
  double *v[3];
  const int nnodes = (int)total_nodes_patch(patch);
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  const int v_ind0 = Ind(subg->comp[0]);
  const int v_ind1 = Ind(subg->comp[1]);
  const int v_ind2 = Ind(subg->comp[2]);
  char *varnames[] = {subg->comp[0],subg->comp[1],subg->comp[2]};
  char desc[MAX_STR_LEN];
  int DB_ret;
  unsigned i,j,k;
  int l;
  
  if (v_ind0 < 0)
    abortEr_s("There is no such field \"%s\" among the fields!\n",subg->comp[0]);
  if (v_ind1 < 0)
    abortEr_s("There is no such field \"%s\" among the fields!\n",subg->comp[1]);
  if (v_ind2 < 0)
    abortEr_s("There is no such field \"%s\" among the fields!\n",subg->comp[2]);
  
  comp[0] = malloc((long unsigned)nnodes*sizeof(*comp[0]));
  pointerEr(comp[0]);
  comp[1] = malloc((long unsigned)nnodes*sizeof(*comp[1]));
  pointerEr(comp[1]);
  comp[2] = malloc((long unsigned)nnodes*sizeof(*comp[2]));
  pointerEr(comp[2]);
  
  v[0] = patch->pool[v_ind0]->v;
  v[1] = patch->pool[v_ind1]->v;
  v[2] = patch->pool[v_ind2]->v;
  
  for (l = 0 ; l < nnodes; ++l)
  {
    IJK((unsigned)l,patch->n,&i,&j,&k);
    unsigned map = i+(unsigned)dims[0]*(j+(unsigned)dims[1]*k);
    comp[0][map] = (float)v[0][l];
    comp[1][map] = (float)v[1][l];
    comp[2][map] = (float)v[2][l];
  }
  
  sprintf(desc,"Vector_%s_%s_%s",
    subg->comp[0],subg->comp[1],subg->comp[2]);
  
  DB_ret = DBPutQuadvar(dbfile,desc,pr->patch->name,3,
    varnames,comp,dims,ndims,0,0,DB_FLOAT,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
    
  /* if there is another file that the field needs to be printed */
  if (pr->file2)
    DB_ret = DBPutQuadvar(pr->file2,desc,pr->patch->name,3,
      varnames,comp,dims,ndims,0,0,DB_FLOAT,DB_NODECENT,0);
    
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
  free(comp[0]);
  free(comp[1]);
  free(comp[2]);
}

/* printing a 3d structured mesh using silo library */
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  DBfile *const dbfile = pr->file;
  float *coords[] = {pr->a,pr->b,pr->c};
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  int DB_ret;
  
  DB_ret = DBPutQuadmesh(dbfile,pr->patch->name,0,coords,dims,ndims,
      DB_FLOAT,DB_NONCOLLINEAR,pr->opt_patch);
      
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
}

/* filling x,y,z value for each node for a structured mesh,
// in compliance with what silo library needs. in fact, silo needs
// x[k][j][i] for each coords and this 3D format maped 
// like i+n0*(j+n1*k) to a 1D format.
*/
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,float **const x,float **const y,float **const z)
{
  float *a, *b, *c;
  const unsigned n0 = patch->n[0];
  const unsigned n1 = patch->n[1];
  const unsigned nn = patch->nn;
  unsigned i,j,k,n;
  
  a = malloc(nn*sizeof(*a));
  pointerEr(a);
  b = malloc(nn*sizeof(*b));
  pointerEr(b);
  c = malloc(nn*sizeof(*c));
  pointerEr(c);
  
  if (strcmp_i(type,"Cartesian"))
  {
    for (n = 0; n < nn; ++n)
    {
      IJK(n,patch->n,&i,&j,&k);
      a[i+n0*(j+n1*k)] = (float)patch->node[n]->x[0];
      b[i+n0*(j+n1*k)] = (float)patch->node[n]->x[1];
      c[i+n0*(j+n1*k)] = (float)patch->node[n]->x[2];
      
    }
  }
  else if (strcmp_i(type,"Curvilinear"))
  {
    for (n = 0; n < nn; ++n)
    {
      IJK(n,patch->n,&i,&j,&k);
      a[i+n0*(j+n1*k)] = (float)patch->node[n]->X[0];
      b[i+n0*(j+n1*k)] = (float)patch->node[n]->X[1];
      c[i+n0*(j+n1*k)] = (float)patch->node[n]->X[2];

    }
  }
  else
    abortEr_s("There is no such type of coordinates"
      " \"%s\"for printing.\n",type);
  
  *x = a;
  *y = b;
  *z = c;
}

/* freeing nodes for silo library */
static void free_nodes_silo(float *x,float *y,float *z)
{
  if (x) free(x);
  if (y) free(y);
  if (z) free(z);
}

/* it reads the parameter to find out what and how it's going to be print.
// given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on.
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {(field1)vs(a,b,c)|coord}.
*/
static void read_parameter_4d(const char *const par,Pr_Field_T *const pr)
{
  struct Info_S *info_s = 0,*Pinfo;
  char *tok = dup_s(par);
  char *save = 0,*sub_tok = 0;
  unsigned Ninfo = 0;
  
  /* check if parameter has been written correctly */
  if (!check_format_s(tok,"?{?}"))
    abortEr(FORMAT_ER_PAR);
  
  sub_tok = sub_s(tok,DL_OC,DL_CC,&save);
  /* => sub_tok = f1,f2,f3,...,(fx,fy,fz),... */
  sub_tok = tok_s(sub_tok,DL_C,&save);
  /* either => sub_tok = f1 and save = f2,f3,... 
  // or     => sub_tok = (fx,fy,fz) and save = ... 
  */
  if (sub_tok == 0)
    abortEr("No field in parameter file is specified for 4d printing.\n");
    
  while (sub_tok)
  {
    char *savess = 0,*ss = 0,*dump;
    
    info_s = realloc(info_s,(Ninfo+1)*sizeof(*info_s)); 
    pointerEr(info_s);
    Pinfo = &info_s[Ninfo];
    Pinfo->field   = 0;
    Pinfo->comp[0] = 0;
    Pinfo->comp[1] = 0;
    Pinfo->comp[2] = 0;
    
    /* if sub_tok = (fx  save = fy,fz)... */
    if (strchr(sub_tok,DL_OP))
    {
      savess = save;
      save = strchr(save,DL_CP);
      save++;
      if (save[0] == DL_C)
        save++;
        
      ss = tok_s(sub_tok,DL_OP,&dump);/* => ss = fx,dump = '\0' */
      Pinfo->comp[0] = dup_s(ss);
      ss = tok_s(0,DL_C,&savess); /* => ss = fy, savess = fz) */
      Pinfo->comp[1] = dup_s(ss);
      ss = tok_s(0,DL_CP,&dump); /* => ss = fz */
      Pinfo->comp[2] = dup_s(ss);
      Pinfo->vec_flg = 1;
    }
    /* sub_tok = f1 */
    else
    {
      Pinfo->field = dup_s(sub_tok);
      Pinfo->vec_flg = 0;
    }
    sub_tok = tok_s(0,DL_C,&save);
    Ninfo++;
  }
  
  pr->group = info_s;
  pr->ng = Ninfo;
  
  free(tok);
}

/* initiating a Pr_Field_T for printing with given grid.
// ->return value: Pr_Field_T
*/
Pr_Field_T *init_PrField(const Grid_T *const grid)
{
  Pr_Field_T *pr = calloc(1,sizeof(*pr));
  pointerEr(pr);
  pr->grid = grid;
  
  return pr;
}

/* freeing Pr_Field_T */
void free_PrField(Pr_Field_T *pr)
{
  free_info_s(pr);
  free(pr);
}
