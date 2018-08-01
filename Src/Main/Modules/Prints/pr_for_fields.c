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
#define DL_B '|'

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
    free(info[i].field);
    free(info[i].axis[0]);
    free(info[i].axis[1]);
    free(info[i].axis[2]);
    free(info[i].coord);
    
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
    DBfile *dbfile = 0;
    char file_name[MAX_STR_LEN];
    float *x = 0,*y = 0,*z = 0;
    float *a = 0,*b = 0,*c = 0;
    unsigned i_pr;
    Flag_T flg;
    
    /* printing Cartesian type */
    flg = NONE;
    for (i_pr = 0; i_pr < npr; ++i_pr)
    {
      if (strcmp_i(pr_info[i_pr].coord,"Cartesian"))
      {
        if (flg == NONE)
        {
          prepare_node_structured_mesh_3d_silo("Cartesian",patch,&x,&y,&z);
          flg = FOUND;
        }
        
        /* opening a file to write in */
        sprintf(file_name,"%s/%s.xyz.%s.%04d.silo",
          pr->folder,pr_info[i_pr].field,patch->name,pr->cycle);
        dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
          "Printing Field on 3D mesh in HDF5 format using silo library",
          DB_HDF5);
        pointerEr(dbfile);
        
        /* setting up some options */
        DBoptlist *opt = DBMakeOptlist(5);
        DBAddOption(opt,DBOPT_XLABEL,(void *)pr_info[i_pr].axis[0]);
        DBAddOption(opt,DBOPT_YLABEL,(void *)pr_info[i_pr].axis[1]);
        DBAddOption(opt,DBOPT_ZLABEL,(void *)pr_info[i_pr].axis[2]);
        DBAddOption(opt,DBOPT_DTIME,&pr->time);
        DBAddOption(opt,DBOPT_CYCLE,&pr->cycle);

        
        /* setting up pr */
        pr->a = x;
        pr->b = y;
        pr->c = z;
        pr->patch = patch;
        pr->opt_patch = opt;
        pr->file = dbfile;
        pr->vptr = &pr_info[i_pr];
        /* printing mesh */
        pr_structured_mesh_3d_silo(pr);
        /* printing filed on mesh */
        pr_field_on_structured_mesh_3d_silo(pr);
        
        DBClose(dbfile);
        DBFreeOptlist(opt);
      }/* end of if (strcmp_i(pr_info[i_pr].coord,"Cartesian")) */
    }/* end of for (i_pr = 0; i_pr < npr; ++i_pr) */
    if (flg == FOUND)
      free_nodes_silo(x,y,z);
    
    /* printing Curvilinear type */
    flg = NONE;
    for (i_pr = 0; i_pr < npr; ++i_pr)
    {
      if (strcmp_i(pr_info[i_pr].coord,"Curvilinear"))
      {
        if (flg == NONE)
        {
          prepare_node_structured_mesh_3d_silo("Curvilinear",patch,&a,&b,&c);
          flg = FOUND;
        }
        
        /* opening a file to write in */
        sprintf(file_name,"%s/%s.abc.%s.%04d.silo",
          pr->folder,pr_info[i_pr].field,patch->name,pr->cycle);
        dbfile = DBCreate(file_name,DB_CLOBBER,DB_LOCAL,
          "Printing Field on 3D mesh in HDF5 format using silo library",
          DB_HDF5);
        pointerEr(dbfile);
        
        /* setting up some options */
        DBoptlist *opt = DBMakeOptlist(5);
        DBAddOption(opt,DBOPT_XLABEL,(void *)pr_info[i_pr].axis[0]);
        DBAddOption(opt,DBOPT_YLABEL,(void *)pr_info[i_pr].axis[1]);
        DBAddOption(opt,DBOPT_ZLABEL,(void *)pr_info[i_pr].axis[2]);
        DBAddOption(opt,DBOPT_DTIME,&pr->time);
        DBAddOption(opt,DBOPT_CYCLE,&pr->cycle);
        
        /* setting up pr */
        pr->a = a;
        pr->b = b;
        pr->c = c;
        pr->patch = patch;
        pr->opt_patch = opt;
        pr->file = dbfile;
        pr->vptr = &pr_info[i_pr];
        /* printing mesh */
        pr_structured_mesh_3d_silo(pr);
        /* printing filed on mesh */
        pr_field_on_structured_mesh_3d_silo(pr);
        
        DBClose(dbfile);
        DBFreeOptlist(opt);
      }/* end of if (strcmp_i(pr_info[i_pr].coord,"Curvilinear")) */
    }/* end of for (i_pr = 0; i_pr < npr; ++i_pr) */
    if (flg == FOUND)
      free_nodes_silo(a,b,c);
    
  }/* end of FOR_ALL_PATCHES(pa,pr->grid) */
}

/* printing field on structured 3d mesh node centered format.
// note data for silo library must be written in cloumn format,
// i.e. 3d to 1d map is i+n0*(j+n1*k).
*/
static void pr_field_on_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
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
   
  DBPutQuadvar1(dbfile,subg->field,pr->patch->name,
    data,dims,ndims,0,0,DB_FLOAT,DB_NODECENT,0);
  
  free(data);
}

/* printing a 3d structured mesh using silo library */
static void pr_structured_mesh_3d_silo(const Pr_Field_T *const pr)
{
  DBfile *dbfile = pr->file;
  float *coords[] = {pr->a,pr->b,pr->c};
  int dims[] = 
    {(int)pr->patch->n[0],(int)pr->patch->n[1],(int)pr->patch->n[2]};
  const int ndims = 3;
  
  DBPutQuadmesh(dbfile,pr->patch->name,0,coords,dims,ndims,
      DB_FLOAT,DB_NONCOLLINEAR,pr->opt_patch);
  
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
  if (!check_format_s(tok,"?{(?)?(?,?,?)|?}?"))
    abortEr(FORMAT_ER_PAR);
  
  
  sub_tok = sub_s(tok,DL_OC,DL_CC,&save);
  /* => sub_tok = (field1)vs(x,y,z)|coord */
  while (sub_tok)
  {
    if (!check_format_s(sub_tok,"(?)?(?,?,?)|?"))
      abortEr(FORMAT_ER_PAR);
  
    char *savef = 0,*sub_tokf = 0;
    char *savec = 0,*sub_tokc = 0;
    char *savess = 0,*ss = 0;
    
    info_s = realloc(info_s,(Ninfo+1)*sizeof(*info_s)); 
    pointerEr(info_s);
    Pinfo = &info_s[Ninfo];
    Pinfo->field = 0;
    
    /* getting field names */
    sub_tokf = sub_s(sub_tok,DL_OP,DL_CP,&savef);/* => sub_tokf = field */
    Pinfo->field = dup_s(sub_tokf);
    
    /* getting axis */
    sub_tokc = sub_s(savef,DL_OP,DL_CP,&savec);/* => sub_tokc = a,b,c */
    ss = tok_s(sub_tokc,DL_C,&savess);/* => ss = a */
    Pinfo->axis[0] = dup_s(ss);
    ss = tok_s(0,DL_C,&savess); /* => ss = b */
    Pinfo->axis[1] = dup_s(ss);
    ss = tok_s(0,DL_C,&savess); /* => ss = c */
    Pinfo->axis[2] = dup_s(ss);
    ss = tok_s(0,DL_B,&savec);
    Pinfo->coord = dup_s(ss);/* ss = coord */
    sub_tok = sub_s(0,DL_OC,DL_CC,&save);
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
