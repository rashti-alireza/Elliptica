/*
// Alireza Rashti
// July 2018
*/

#include "pr_for_fields_hdf5.h"


/* Print Fields:
// note: in Silo language, I'm using curvilinear format for both mesh and data (fields)
// regardless of the patch is Cartesian or not.
// note: data and mesh in Silo must be written in column major order,
//       otherwise for inhomogeneous resolutions you get a scrambled mesh.
//
// usage examples:
// ===============
// # parameter that is determined in input file is like (supporting regular expression):
// silo_output_3d = (V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2),^dpsi_D.+
// # as one can see the vector quantities determined by parenthesis 
//
// Pr_Field_T *pr  = init_PrField(grid);
// pr->folder      = "folder_path";
// pr->cycle       = iteration_number;// if you wanna plot data at each iteration
//
// # the following options and flag are not necessary, their default value is 0.
// pr->multimesh_f = 1; # if you wanna make a master file for all patches as a whole grid
// pr->multivar_f  = 1; # if you wanna make a master file for all fields. This, option seems not working with VisIt.
// pr->abc_f       = 1; # if you wanna have the patches and fields in (X,Y,Z) corods or (a,b,c) coords.
//
// # print all patchs and fields
// pr_fields(pr);
//
// # free
// free_PrField(pr);
*/

/* row major to column major order.
// NOTE: silo using column major format thus if the fields saved in 
// row major format, which is the default, one must convert them. */
#define row2col(i,j,k) (i+n[0]*(j+n[1]*k))

/* printing fields with HDF5 format using silo library.
// note: here patch and mesh are used interchangeably and grid
// composed of all patches (meshes). */
void pr_hdf5_silo(Pr_Field_T *const pr)
{
  struct Info_S *const pr_info = pr->group;
  const Uint npr = pr->ng;
  Uint pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    DBfile *dbfile_xyz = 0;/* file for cartesian  */
    Uint i_pr;
   
    /* printing in Cartesian coords, namely in (x,y,z) values */
    dbfile_xyz = make_structured_mesh_3d_xyz(pr,patch);
    
    /* printing field on the meshes */
    for (i_pr = 0; i_pr < npr; ++i_pr)
    {
      pr->file  = dbfile_xyz;
      pr->vptr  = &pr_info[i_pr];
      
      /* if it is vector to be printed */
      if (pr_info[i_pr].vec_flg)
        pr_vector_on_structured_mesh_3d_silo(pr);
      
      /* printing scalar field on mesh */
      else
        pr_scalar_on_structured_mesh_3d_silo(pr);
    }
    
    if (DBClose(dbfile_xyz) == -1)
      Error0("Silo library failed to close the file.\n");
  }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
  
  /* if you want in (a,b,c) coords */
  if (pr->abc_f)
  {
    FOR_ALL_PATCHES(pa,pr->grid)
    {
      Patch_T *patch = pr->grid->patch[pa];
      DBfile *dbfile_abc = 0;/* file for curvilinear */
      Uint i_pr;
     
      /* printing Curvilinear coords, namely in (a,b,c) values */
      dbfile_abc = make_structured_mesh_3d_abc(pr,patch);
      
      /* printing field on the meshes */
      for (i_pr = 0; i_pr < npr; ++i_pr)
      {
        pr->file  = dbfile_abc;
        pr->vptr  = &pr_info[i_pr];
        
        /* if it is vector to be printed */
        if (pr_info[i_pr].vec_flg)
          pr_vector_on_structured_mesh_3d_silo(pr);
        
        /* printing scalar field on mesh */
        else
          pr_scalar_on_structured_mesh_3d_silo(pr);
      }
      
      if (DBClose(dbfile_abc) == -1)
        Error0("Silo library failed to close the file.\n");
    }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
  
  }
    
  if (pr->multimesh_f)/* if multimesh flag active */
    write_multi_mesh(pr);
  if (pr->multivar_f)/* if multivar flag active */
    write_multi_vars(pr);
}

/* make master file to encompass all of the patch info 
// for convenience of visulization, so one doesn't need to 
// add each patch (mesh) one by one. */
static void write_multi_mesh(const Pr_Field_T *const pr)
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
  Uint pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    
    /* populating DBPutMultimesh function arguments */
    patch_names_xyz[pa] = calloc(MAX_STR_LEN,1);
    IsNull(patch_names_xyz[pa]);
    sprintf(patch_names_xyz[pa],"%s_xyz_%04d.silo:%s",/* this is how we set patch files */
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
  IsNull(grid_file_xyz);
  sprintf(grid_name,"grid%d_xyz",pr->grid->gn);
  DB_ret = DBPutMultimesh
    (grid_file_xyz,grid_name,npatch,
    (char const * const *)patch_names_xyz,patch_types,0);
  if (DB_ret == -1)
    Error0("Silo library failed to make multi-mesh.\n");
    
  /* close and free */
  if (DBClose(grid_file_xyz) == -1)
    Error0("Silo library failed to close the file.\n");
  for (i = 0; i < npatch; ++i)
    free(patch_names_xyz[i]);
  
  /* if you want in (a,b,c) coords */
  if (pr->abc_f)
  {
    FOR_ALL_PATCHES(pa,pr->grid)
    {
      Patch_T *patch = pr->grid->patch[pa];
      
      /* populating DBPutMultimesh function arguments */
      patch_names_abc[pa] = calloc(MAX_STR_LEN,1);
      IsNull(patch_names_abc[pa]);
      sprintf(patch_names_abc[pa],"%s_abc_%04d.silo:%s",/* this is how we set patch files */
                               patch->name,pr->cycle,patch->name);
      patch_types[pa] = DB_QUAD_CURV;/* note, in my implementation both meshes 
                                     // are curvilinear type which makes it more flexible */
    }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
    
    /* make multi-mesh master file for abc type */
    sprintf(grid_file_path_abc,"%s/grid%d_abc_%04d.silo",
                           pr->folder,pr->grid->gn,pr->cycle);
    grid_file_abc = 
      DBCreate(grid_file_path_abc,DB_CLOBBER,DB_LOCAL,
              "Master file composed of all patches path",DB_HDF5);
    IsNull(grid_file_abc);
    sprintf(grid_name,"grid%d_abc",pr->grid->gn);
    DB_ret = DBPutMultimesh
      (grid_file_abc,grid_name,npatch,
      (char const * const *)patch_names_abc,patch_types,0);
    if (DB_ret == -1)
      Error0("Silo library failed to make multi-mesh.\n");
      
    /* close and free */
    if (DBClose(grid_file_abc) == -1)
      Error0("Silo library failed to close the file.\n");
    for (i = 0; i < npatch; ++i)
      free(patch_names_abc[i]);
  }/* end of if (pr->abc_f) */
}

/* make master file to encompass ALL of the VARIABLES info 
// for convenience of visulization, so one doesn't need to 
// add each variables one by one. */
static void write_multi_vars(const Pr_Field_T *const pr)
{
  struct Info_S *const pr_info = pr->group;
  const Uint npr = pr->ng;
  Uint i_pr;
  
  /* for all fields or variables */
  for (i_pr = 0; i_pr < npr; ++i_pr)
  {
    /* if it is vector */
    if (pr_info[i_pr].vec_flg)
    {
      make_multi_var(pr,pr_info[i_pr].comp[0]);
      make_multi_var(pr,pr_info[i_pr].comp[1]);
      make_multi_var(pr,pr_info[i_pr].comp[2]);
    }
    /* if it is scalar */
    else
      make_multi_var(pr,pr_info[i_pr].field);
  }

}

/* make master file for ONE variable to encompass all of the variable info
// at all patches for convenience of visulization, so one doesn't need to 
// add each patch one by one. */  
static void make_multi_var(const Pr_Field_T *const pr,const char *const var)
{
  const int npatch = (int)pr->grid->np;
  DBfile *var_file_xyz  = 0;
  DBfile *var_file_abc  = 0;
  char *var_names_xyz[npatch];
  char *var_names_abc[npatch];
  char var_file_path_xyz[MAX_STR_LEN];
  char var_file_path_abc[MAX_STR_LEN];
  char var_name[MAX_STR_LEN];
  int  var_types[npatch];
  int DB_ret,i;
  Uint pa;
  
  FOR_ALL_PATCHES(pa,pr->grid)
  {
    Patch_T *patch = pr->grid->patch[pa];
    
    if (_Ind(var) < 0)
      continue;
    
    /* populating DBPutMultimesh function arguments */
    var_names_xyz[pa] = calloc(MAX_STR_LEN,1);
    IsNull(var_names_xyz[pa]);
    sprintf(var_names_xyz[pa],"%s_xyz_%04d.silo:%s",/* this is how we set var files */
                             patch->name,pr->cycle,var);
    var_types[pa] = DB_QUADVAR;
  }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
  
  /* make multi-var master file for xyz type */
  sprintf(var_file_path_xyz,"%s/grid%d_%s_xyz_%04d.silo",
                            pr->folder,pr->grid->gn,var,pr->cycle);
  var_file_xyz = 
    DBCreate(var_file_path_xyz,DB_CLOBBER,DB_LOCAL,
            "Master file composed of all fields path",DB_HDF5);
  IsNull(var_file_xyz);
  sprintf(var_name,"%s_xyz",var);
  DB_ret = DBPutMultivar
    (var_file_xyz,var_name,npatch,
    (char const * const *)var_names_xyz,var_types,0);
  if (DB_ret == -1)
    Error0("Silo library failed to make multi-var.\n");
    
  /* close and free */
  if (DBClose(var_file_xyz) == -1)
    Error0("Silo library failed to close the file.\n");
  for (i = 0; i < npatch; ++i)
    free(var_names_xyz[i]);

//test
  //IsNull(DBGetMultivar(var_file_xyz,var_name));
//end

  /* if you want in (a,b,c) coords */
  if (pr->abc_f)
  {
    FOR_ALL_PATCHES(pa,pr->grid)
    {
      Patch_T *patch = pr->grid->patch[pa];
      
      if (_Ind(var) < 0)
        continue;
      
      /* populating DBPutMultimesh function arguments */
      var_names_abc[pa] = calloc(MAX_STR_LEN,1);
      IsNull(var_names_abc[pa]);
      sprintf(var_names_abc[pa],"%s_abc_%04d.silo:%s",/* this is how we set var files */
                               patch->name,pr->cycle,var);
      var_types[pa] = DB_QUADVAR;
    }/* end of FOR_ALL_PATCHES(pa,pr->grid) */

    /* make multi-var master file for abc type */
    sprintf(var_file_path_abc,"%s/grid%d_%s_abc_%04d.silo",
                              pr->folder,pr->grid->gn,var,pr->cycle);
    var_file_abc = 
      DBCreate(var_file_path_abc,DB_CLOBBER,DB_LOCAL,
              "Master file composed of all fields path",DB_HDF5);
    IsNull(var_file_abc);
    sprintf(var_name,"%s_abc",var);
    DB_ret = DBPutMultivar
      (var_file_abc,var_name,npatch,
      (char const * const *)var_names_abc,var_types,0);
    if (DB_ret == -1)
      Error0("Silo library failed to make multivar.\n");
      
    /* close and free */
    if (DBClose(var_file_abc) == -1)
      Error0("Silo library failed to close the file.\n");
    for (i = 0; i < npatch; ++i)
      free(var_names_abc[i]);
  }
}

/* printing 3D mesh in Cartesian coordinates
// ->return value: the file containing the mesh. */
static void *make_structured_mesh_3d_xyz(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  const Uint nn = patch->nn;
  double *x = alloc_double(nn),
         *y = alloc_double(nn),
         *z = alloc_double(nn);
  char *label[3];
  const char *mesh_name = strstr(patch->name,"_");/* grid\d?_ */
  assert(mesh_name);
  mesh_name++;
  
  prepare_node_structured_mesh_3d_silo("Cartesian",patch,x,y,z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s_xyz_%04d.silo",
    pr->folder,mesh_name,pr->cycle);
  dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
    "3D mesh in Cartesian values with HDF5 format using silo library",
    DB_HDF5);
  IsNull(dbfile);

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
static void *make_structured_mesh_3d_abc(Pr_Field_T *const pr,const Patch_T *const patch)
{
  DBfile *dbfile = 0;
  char file_name[MAX_STR_LEN];
  const Uint nn = patch->nn;
  double *x = alloc_double(nn),
         *y = alloc_double(nn),
         *z = alloc_double(nn);
  char *label[3];
  const char *mesh_name = strstr(patch->name,"_");/* grid\d?_ */
  assert(mesh_name);
  mesh_name++;
  
  prepare_node_structured_mesh_3d_silo("Curvilinear",patch,x,y,z);
  
  /* opening a file to write in for Cartesian */
  sprintf(file_name,"%s/%s_abc_%04d.silo",
    pr->folder,mesh_name,pr->cycle);
  dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
    "3D Curvilinear mesh with HDF5 format using silo library",
    DB_HDF5);
  IsNull(dbfile);

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
  const Uint nn = patch->nn;
  const Uint *const n = patch->n;
  struct Info_S *subg = pr->vptr;
  double *data = 0;
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  int v_ind = _Ind(subg->field);
  int DB_ret;
  const char *mesh_name = strstr(pr->patch->name,"_");/* grid\d?_ */
  Uint ijk,i,j,k;
  
  assert(mesh_name);
  mesh_name++;
  
  /* if field did not find a match maybe there is a regex */
  if (v_ind < 0)
  {
    Uint Nf   = 0;
    Uint *fInd = find_field_index_regex(patch,subg->field,&Nf);
    Uint f = 0;
    
    for (f = 0; f < Nf; ++f)
    {
      if (!patch->fields[fInd[f]]->v)
        continue;
        
      /* fields value */
      data = alloc_double(nn);
      for (ijk = 0; ijk < nn; ++ijk)
      {
        ijk_to_i_j_k(ijk,n,&i,&j,&k);
        data[row2col(i,j,k)] = patch->fields[fInd[f]]->v[ijk];
      }
       
      DB_ret = DBPutQuadvar1(dbfile,patch->fields[fInd[f]]->name,mesh_name,
        data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
      if (DB_ret == -1)
        Error0("Silo library failed to print.\n");
      
      /* if there is another file that the field needs to be printed */
      if (pr->file2)
      {
        DB_ret = DBPutQuadvar1(pr->file2,patch->fields[fInd[f]]->name,mesh_name,
          data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
        if (DB_ret == -1)
          Error0("Silo library failed to print.\n");
      }
      
      Free(data);
    }
    Free(fInd);
  }
  else
  {
    if (!patch->fields[v_ind]->v)
      return;
      
    /* fields value */
    data = alloc_double(nn);
    for (ijk = 0; ijk < nn; ++ijk)
    {
      ijk_to_i_j_k(ijk,n,&i,&j,&k);
      data[row2col(i,j,k)] = patch->fields[v_ind]->v[ijk];
    }
     
    DB_ret = DBPutQuadvar1(dbfile,subg->field,mesh_name,
      data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
    if (DB_ret == -1)
      Error0("Silo library failed to print.\n");
    
    /* if there is another file that the field needs to be printed */
    if (pr->file2)
    {
      DB_ret = DBPutQuadvar1(pr->file2,subg->field,mesh_name,
        data,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
      if (DB_ret == -1)
        Error0("Silo library failed to print.\n");
    }
    
    Free(data);
  }/* else */
}

/* printing vector field on structured 3d mesh with node centered format. */
static void pr_vector_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  assert(pr->file);
  
  DBfile *const dbfile = pr->file; 
  const Patch_T *const patch = pr->patch;
  const Uint nn = patch->nn;
  const Uint *const n = patch->n;
  struct Info_S *subg = pr->vptr;
  double *comp[3] = {0};
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  const int v_ind0 = _Ind(subg->comp[0]);
  const int v_ind1 = _Ind(subg->comp[1]);
  const int v_ind2 = _Ind(subg->comp[2]);
  char *varnames[] = {subg->comp[0],subg->comp[1],subg->comp[2]};
  char desc[MAX_STR_LEN];
  int DB_ret;
  const char *mesh_name = strstr(pr->patch->name,"_");/* grid\d?_ */
  Uint ijk,i,j,k;
  
  assert(mesh_name);
  mesh_name++;
  
  /* if field does not exists or empty */
  if (
      v_ind0 < 0 || 
      v_ind1 < 0 || 
      v_ind2 < 0 ||
      !patch->fields[v_ind0]->v ||
      !patch->fields[v_ind1]->v || 
      !patch->fields[v_ind2]->v
      )
    return;
  
  comp[0] = alloc_double(nn);
  comp[1] = alloc_double(nn);
  comp[2] = alloc_double(nn);
  
  /* change the order from row major to column major order */
  for (ijk = 0; ijk < nn; ++ijk)
  {
    ijk_to_i_j_k(ijk,n,&i,&j,&k);
    comp[0][row2col(i,j,k)] = patch->fields[v_ind0]->v[ijk];
    comp[1][row2col(i,j,k)] = patch->fields[v_ind1]->v[ijk];
    comp[2][row2col(i,j,k)] = patch->fields[v_ind2]->v[ijk];
  }
  
  sprintf(desc,"Vector_%s_%s_%s",
    subg->comp[0],subg->comp[1],subg->comp[2]);
  
  DB_ret = DBPutQuadvar
    (dbfile,desc,mesh_name,3,(char const * const *)varnames,
    comp,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
  if (DB_ret == -1)
    Error0("Silo library failed to print.\n");
    
  /* if there is another file that the field needs to be printed */
  if (pr->file2)
  {
    DB_ret = DBPutQuadvar
      (pr->file2,desc,mesh_name,3,(char const * const *)varnames,
       comp,dims,ndims,0,0,DB_DOUBLE,DB_NODECENT,0);
    if (DB_ret == -1)
      Error0("Silo library failed to print.\n");
  }
  
  Free(comp[0]);
  Free(comp[1]);
  Free(comp[2]);
  
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
  const char *mesh_name = strstr(pr->patch->name,"_");/* grid\d?_ */
  assert(mesh_name);
  mesh_name++;
  
  DB_ret = DBPutQuadmesh(dbfile,mesh_name,0,coords,dims,ndims,
      DB_DOUBLE,DB_NONCOLLINEAR,pr->opt_patch);
      
  if (DB_ret == -1)
    Error0("Silo library failed to print.\n"); 
}

/* filling x,y,z value for each node for a structured mesh,
// in compliance with what silo library needs. in fact, silo needs
// x[i][j][k] for each coords. */
static void prepare_node_structured_mesh_3d_silo(const char *const type,const Patch_T *const patch,double *const x,double *const y,double *const z)
{
  const Uint nn = patch->nn;
  const Uint *const n = patch->n;
  Uint ijk,i,j,k;
  
  if (strcmp_i(type,"Cartesian"))
  {
    for (ijk = 0; ijk < nn; ++ijk)
    {
      ijk_to_i_j_k(ijk,patch->n,&i,&j,&k);
      
      x[row2col(i,j,k)] = patch->node[ijk]->x[0];
      y[row2col(i,j,k)] = patch->node[ijk]->x[1];
      z[row2col(i,j,k)] = patch->node[ijk]->x[2];
    }
  }
  else if (strcmp_i(type,"Curvilinear"))
  {
    for (ijk = 0; ijk < nn; ++ijk)
    {
      ijk_to_i_j_k(ijk,patch->n,&i,&j,&k);
      
      x[row2col(i,j,k)] = patch->node[ijk]->X[0];
      y[row2col(i,j,k)] = patch->node[ijk]->X[1];
      z[row2col(i,j,k)] = patch->node[ijk]->X[2];
    }
  }
  else
    Errors("There is no such type of coordinates"
      " \"%s\"for printing.\n",type);
  
}
