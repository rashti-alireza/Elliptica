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
// # parameter that is determined in input file is like:
// print_fields_4d = yes,Format:HDF5,{(V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2)}
// # as one can see the vector quantities determined by parenthesis 
// # and all of the desired fields need to be put in curly bracket
//
// Pr_Field_T *pr  = init_PrField(grid);
// pr->folder      = "folder_path";
// pr->par         = "print_fields_4d";
// pr_fields(pr);
// free_PrField(pr);
*/

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and prints the whole grid and fields in the specified folder. */
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

/* printing fields with HDF5 format using silo library.
// note: here patch and mesh are used interchangeably and grid
// composed of all patches (meshes). */
static void pr_fields_on_grid_HDF5_4d(Pr_Field_T *const pr)
{
  struct Info_S *const pr_info = pr->group;
  const unsigned npr = pr->ng;
  unsigned pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    DBfile *dbfile_Cart = 0;/* file cartesian  */
    DBfile *dbfile_Curv = 0;/* file for curvilinear */
    unsigned i_pr;
   
    /* printing Cartesian mesh, 
    // namely coords are Cartesian, e.g (x,y,z) */
    dbfile_Cart = make_structured_mesh_3d_Cartesian(pr,patch);
    /* printing Curvilinear mesh, 
    // namely coords are Curvilinear, e.g (r,theta,phi) */
    dbfile_Curv = make_structured_mesh_3d_Curvilinear(pr,patch);
    
    /* printing field on the meshes */
    for (i_pr = 0; i_pr < npr; ++i_pr)
    {
      pr->file  = dbfile_Cart;
      pr->file2 = dbfile_Curv;
      pr->vptr  = &pr_info[i_pr];
      
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
    
  pr_multi_mesh(pr);
}

/* make master file to encompass all of the patch info 
// for convenience of visulization, so one doesn't need to add each box separately */
static void pr_multi_mesh(Pr_Field_T *const pr)
{
  const int npatch   = (int)pr->grid->np;
  DBfile *grid_file_xyz  = 0;
  DBfile *grid_file_abc  = 0;
  char *patch_names_xyz[npatch];
  char *patch_names_abc[npatch];
  char  grid_file_path_xyz[MAX_STR_LEN];/* for Cartesian coords  */
  char  grid_file_path_abc[MAX_STR_LEN];/* for Curvilinear coords */
  char  grid_name[MAX_STR_LEN];
  int   patch_types[npatch];
  int DB_ret,i; 
  unsigned pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    
    /* populating DBPutMultimesh function arguments */
    patch_names_xyz[pa] = calloc(MAX_STR_LEN,1);
    pointerEr(patch_names_xyz[pa]);
    patch_names_abc[pa] = calloc(MAX_STR_LEN,1);
    pointerEr(patch_names_abc[pa]);
    sprintf(patch_names_xyz[pa],"%s_xyz_%04d.silo:%s",/* this is how we set patch files */
                             patch->name,pr->cycle,patch->name);
    sprintf(patch_names_abc[pa],"%s_abc_%04d.silo:%s",/* this is how we set patch files */
                             patch->name,pr->cycle,patch->name);
    patch_types[pa] = DB_QUAD_CURV;/* note, in my implementation both meshes 
                                   // are curvilinear type which makes it more flexible */
  }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
  
  /* make multi-mesh master file for xyz type */
  sprintf(grid_file_path_xyz,"%s/grid%d_xyz_%04d.silo",
                         pr->folder,pr->grid->gn,pr->cycle);
  grid_file_xyz = 
    DBCreate(grid_file_path_xyz,DB_CLOBBER,DB_LOCAL,
            "Master file composed of all patches path",DB_HDF5);
  pointerEr(grid_file_xyz);
  sprintf(grid_name,"grid%d_xyz",pr->grid->gn);
  DB_ret = DBPutMultimesh(grid_file_xyz,grid_name,npatch,patch_names_xyz,patch_types,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to make multi-mesh.\n");
    
  /* close and free */
  if (DBClose(grid_file_xyz) == -1)
    abortEr("Silo library failed to close the file.\n");
  for (i = 0; i < npatch; ++i)
    free(patch_names_xyz[i]);

  /* make multi-mesh master file for abc type */
  sprintf(grid_file_path_abc,"%s/grid%d_abc_%04d.silo",
                         pr->folder,pr->grid->gn,pr->cycle);
  grid_file_abc = 
    DBCreate(grid_file_path_abc,DB_CLOBBER,DB_LOCAL,
            "Master file composed of all patches path",DB_HDF5);
  pointerEr(grid_file_abc);
  sprintf(grid_name,"grid%d_abc",pr->grid->gn);
  DB_ret = DBPutMultimesh(grid_file_abc,grid_name,npatch,patch_names_abc,patch_types,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to make multi-mesh.\n");
    
  /* close and free */
  if (DBClose(grid_file_abc) == -1)
    abortEr("Silo library failed to close the file.\n");
  for (i = 0; i < npatch; ++i)
    free(patch_names_abc[i]);
    
}

/* printing 3D mesh in Cartesian coordinates
// ->return value: the file containing the mesh. */
static void *make_structured_mesh_3d_Cartesian(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  const unsigned nn = patch->nn;
  double *x = alloc_double(nn),
         *y = alloc_double(nn),
         *z = alloc_double(nn);
  char *label[3];
  prepare_node_structured_mesh_3d_silo("Cartesian",patch,x,y,z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s_xyz_%04d.silo",
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
  
  /* freeing */
  DBFreeOptlist(opt);
  free(x);
  free(y);
  free(z);
  free(label[0]);
  free(label[1]);
  free(label[2]);
  
  return dbfile;
}

/* printing 3D in Curvilinear coordinates
// ->return value: the file containing the mesh.
*/
static void *make_structured_mesh_3d_Curvilinear(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  const unsigned nn = patch->nn;
  double *x = alloc_double(nn),
         *y = alloc_double(nn),
         *z = alloc_double(nn);
  char *label[3];
  prepare_node_structured_mesh_3d_silo("Curvilinear",patch,x,y,z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s_abc_%04d.silo",
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
  
  /* freeing */
  DBFreeOptlist(opt);
  free(x);
  free(y);
  free(z);
  free(label[0]);
  free(label[1]);
  free(label[2]);

  return dbfile;
}

/* printing scalar on structured 3d mesh with node centered format. */
static void pr_scalar_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  assert(pr->file);
  
  DBfile *const dbfile = pr->file; 
  const Patch_T *const patch = pr->patch;
  struct Info_S *subg = pr->vptr;
  double *data;
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  const int v_ind = _Ind(subg->field);
  int DB_ret;
  
  if (v_ind < 0)
    return;
  
  /* fields value */
  data = patch->pool[v_ind]->v;
   
  DB_ret = DBPutQuadvar1(dbfile,subg->field,pr->patch->name,
    data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
  /* if there is another file that the field needs to be printed */
  if (pr->file2)
    DB_ret = DBPutQuadvar1(pr->file2,subg->field,pr->patch->name,
      data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
}

/* printing vector field on structured 3d mesh with node centered format. */
static void pr_vector_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  assert(pr->file);
  
  DBfile *const dbfile = pr->file; 
  const Patch_T *const patch = pr->patch;
  struct Info_S *subg = pr->vptr;
  double *comp[3];
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  const int v_ind0 = _Ind(subg->comp[0]);
  const int v_ind1 = _Ind(subg->comp[1]);
  const int v_ind2 = _Ind(subg->comp[2]);
  char *varnames[] = {subg->comp[0],subg->comp[1],subg->comp[2]};
  char desc[MAX_STR_LEN];
  int DB_ret;
  
  if (v_ind0 < 0 || v_ind1 < 0 || v_ind2 < 0)
    return;
  
  comp[0] = patch->pool[v_ind0]->v;
  comp[1] = patch->pool[v_ind1]->v;
  comp[2] = patch->pool[v_ind2]->v;
  
  sprintf(desc,"Vector_%s_%s_%s",
    subg->comp[0],subg->comp[1],subg->comp[2]);
  
  DB_ret = DBPutQuadvar(dbfile,desc,pr->patch->name,3,
    varnames,comp,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
    
  /* if there is another file that the field needs to be printed */
  if (pr->file2)
    DB_ret = DBPutQuadvar(pr->file2,desc,pr->patch->name,3,
      varnames,comp,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
    
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n");
  
}

/* printing a 3d structured mesh using silo library */
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  DBfile *const dbfile = pr->file;
  const double *coords[] = {pr->a,pr->b,pr->c};
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  int DB_ret;
  
  DB_ret = DBPutQuadmesh(dbfile,pr->patch->name,0,coords,dims,ndims,
      DB_DOUBLE,DB_NONCOLLINEAR,pr->opt_patch);
      
  if (DB_ret == -1)
    abortEr("Silo library failed to print.\n"); 
}

/* filling x,y,z value for each node for a structured mesh,
// in compliance with what silo library needs. in fact, silo needs
// x[i][j][k] for each coords. */
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,double *const x,double *const y,double *const z)
{
  const unsigned nn = patch->nn;
  unsigned ijk;
  
  if (strcmp_i(type,"Cartesian"))
  {
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[ijk] = patch->node[ijk]->x[0];
      y[ijk] = patch->node[ijk]->x[1];
      z[ijk] = patch->node[ijk]->x[2];
    }
  }
  else if (strcmp_i(type,"Curvilinear"))
  {
    for (ijk = 0; ijk < nn; ++ijk)
    {
      x[ijk] = patch->node[ijk]->X[0];
      y[ijk] = patch->node[ijk]->X[1];
      z[ijk] = patch->node[ijk]->X[2];
    }
  }
  else
    abortEr_s("There is no such type of coordinates"
      " \"%s\"for printing.\n",type);
  
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
