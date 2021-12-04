/*
// Alireza Rashti
// June 2018
*/

#include "pr_for_tests.h"
#define MAXSTR 400
/* printing different quantities for test */
/* ->return value: if print option is on return 1 otherwise 0 */
int test_print(const Print_T f)
{
  const char *on; 
  
  switch(f)
  {
    case PRINT_PARAMETERS:
      on = PgetsEZ("print_parameters");
      if (on == 0) return 0;
      if (strcmp_i(on,"yes")|| strcmp_i(on,"y"))
        return 1;
      break;
    case PRINT_COORDS:
      on = PgetsEZ("print_coords");
      if (on == 0) return 0;
      if (strcmp_i(on,"yes")|| strcmp_i(on,"y"))
        return 1;
      break;
    case PRINT_INTERFACES:
      on = PgetsEZ("print_interfaces");
      if (on == 0) return 0;
      if (strcmp_i(on,"yes")|| strcmp_i(on,"y"))
        return 1;
      break;  
    default:
      break;
  }
  
  return 0;
}

/* print interfaces */
void pr_interfaces(const Grid_T *const grid)
{
  FILE *f;
  char str[MAXSTR]={'\0'},*path;
  const char *path_par;
  Interface_T **face;
  Node_T **node;
  SubFace_T *subf,*subf2;
  enum archive_e 
  {
    C = 0/* copy pair */,
    IN_P/* interpolation pair */,
    O/* outerbound */,
    IN/* innerboundary */,
    EX/* external boundary */,
    INT/* internal boundary */,
    IN_UNP/* interpolation unpair */,
    TOT_ARCH/* total number of above */
  };
  struct Archive_S *arch[TOT_ARCH];
  Uint N[TOT_ARCH];
  Uint pa,fc,sf,i,j,n,save;
  Flag_T flg;
  
  /* initializing */
  for (i = 0; i < TOT_ARCH; ++i)
  {
    arch[i] = 0;
    N[i] = 0;
  }
  
  
  if (get_parameter("Diagnostics"))
    path_par = PgetsEZ("Diagnostics");
  else
    path_par = PgetsEZ("top_directory");
  
  path = make_directory(path_par,"InterfaceInfo");
  
  str[0] = '\0';
  sprintf(str,"%s/interface_pairings.rm",path);
  f = Fopen(str,"w");
  
  /* printing out pairing information: outerB*/
  fprintf(f,"###outer-boundaries are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->outerB == 1)
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u(%s),face:%u,subface:%u\n  '%s'\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          
          add_to_archive(&arch[O],subf,0,&N[O],"outerBndry");
          
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
  fprintf(f,"%s\n",PR_LINE);  
  fflush(f);

  /* printing out pairing information: innerB*/
  fprintf(f,"###inner-boundaries are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->innerB == 1)
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n  '%s'\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          
          add_to_archive(&arch[IN],subf,0,&N[IN],"innerBndry");
          
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
  fprintf(f,"%s\n",PR_LINE);
  fflush(f);

  /* printing out pairing information: internal face*/
  fprintf(f,"###internal face are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->exterF == 0)
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n  '%s'\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          
          add_to_archive(&arch[INT],subf,0,&N[INT],"internalBndy");
          
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
  fprintf(f,"%s\n",PR_LINE);
  fflush(f);
  
  /* printing out pairing information: interpolation inside a patch */
  fprintf(f,"###interploation interfaces happen inside a patch are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->exterF == 1 && subf->outerB == 0 && 
            subf->innerB == 0 && subf->touch == 0 )
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n  '%s'\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          
          add_to_archive(&arch[IN_UNP],subf,0,&N[IN_UNP],"innerIntrpltn");
          
          flg = FOUND;
        }
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
  fprintf(f,"%s\n",PR_LINE);
  fflush(f);
    
  /* printing out pairing information: interpolation*/
  fprintf(f,"###interpolation interfaces are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if ( subf->copy == 0 && subf->touch == 1 )
        {
          subf2 = get_paired_subface(subf);
          
          save = N[IN_P];
          add_to_archive(&arch[IN_P],subf,subf2,&N[IN_P],"faceIntrptn");
          
          /* avoiding double writing */
          if (save != N[IN_P])
          {
            str[0] = '\0';
            sprintf(str,"|->faceIntrptn%u:a{ %s }\n",N[IN_P]-1,subf->flags_str);
            fprintf(f,"%s",str);
            str[0] = '\0';
            sprintf(str,"|->faceIntrptn%u:b{ %s }\n",N[IN_P]-1,subf2->flags_str);
            fprintf(f,"%s",str);
            fprintf(f,"%s\n",PR_LINE);
          }
          
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
  fprintf(f,"%s\n",PR_LINE);
  fflush(f);
  
  /* printing out pairing information: copy*/
  fprintf(f,"###copying interfaces are:\n");
  flg = NONE;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->copy == 1)
        {
          subf2 = get_paired_subface(subf);
          
          save = N[C];
          add_to_archive(&arch[C],subf,subf2,&N[C],"faceCpy");
          
          /* avoiding double writing */
          if (save != N[C])
          {
            str[0] = '\0';
            sprintf(str,"|->faceCpy%u:a{ %s }\n",N[C]-1,subf->flags_str);
            fprintf(f,"%s",str);
            str[0] = '\0';
            sprintf(str,"|->faceCpy%u:b{ %s }\n",N[C]-1,subf2->flags_str);
            fprintf(f,"%s",str);
            fprintf(f,"%s\n",PR_LINE);
          }

          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
  
  fprintf(f,"%s\n",PR_LINE);
  /* count Neumann or Dirichlet for each patch
  // note: it's supposed each interface has only one type of B.C. . */
  int Neumann_tot   = 0;
  int Dirichlet_tot = 0;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    int Neumann_count   = 0;
    int Dirichlet_count = 0;
    
    FOR_ALL(fc,face)
    {
      if (face[fc]->df_dn_set)
      {
        if (face[fc]->df_dn)
          Neumann_count++;
        else
          Dirichlet_count++;
      }
    }
    Neumann_tot  += Neumann_count;
    Dirichlet_tot+= Dirichlet_count;
    
    if (1)
    printf("%-43s = {Neumann = %d, Dirichlet = %d}\n",
           grid->patch[pa]->name,Neumann_count,Dirichlet_count);
    fprintf(f,"%-43s = {Neumann = %d, Dirichlet = %d}\n",
           grid->patch[pa]->name,Neumann_count,Dirichlet_count);
  
  }
  if(1)
  printf("total Neumann = %d, total Dirichlet = %d\n",
          Neumann_tot,Dirichlet_tot);
  fprintf(f,"total Neumann = %d, total Dirichlet = %d\n",
          Neumann_tot,Dirichlet_tot);

  Fclose(f);
  
  /* printing interfacess */
  for (n = 0; n < TOT_ARCH; ++n)
  {
    for (i = 0; i < N[n]; ++i)
    {
      if (arch[n][i].s1 && arch[n][i].s2)
      {
        str[0] = '\0';
        sprintf(str,"%s/%s.paired",path,arch[n][i].n1);
        f = Fopen(str,"w");
        fprintf(f,"#%s\n",arch[n][i].s1->flags_str);
        
        node = arch[n][i].s1->patch->node;
        
        for (j = 0; j < arch[n][i].s1->np; ++j)
        {
          double *x = node[arch[n][i].s1->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        Fclose(f);
        
        str[0] = '\0';
        sprintf(str,"%s/%s.paired",path,arch[n][i].n2);
        f = Fopen(str,"w");
        fprintf(f,"#%s\n",arch[n][i].s2->flags_str);
        
        node = arch[n][i].s2->patch->node;
        
        for (j = 0; j < arch[n][i].s2->np; ++j)
        {
          double *x = node[arch[n][i].s2->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        Fclose(f);
      }/* if (arch[n][i].s1 && arch[n][i].s2) */
      else
      {
        str[0] = '\0';
        sprintf(str,"%s/%s.single",path,arch[n][i].n1);
        f = Fopen(str,"w");
        fprintf(f,"#%s\n",arch[n][i].s1->flags_str);
        
        node = arch[n][i].s1->patch->node;
        
        for (j = 0; j < arch[n][i].s1->np; ++j)
        {
          double *x = node[arch[n][i].s1->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        Fclose(f);
      
      }
    }/* end of for (i = 0; i < N[n]; ++i) */
  }/* end of for (j = 0; j < TOT_ARCH; ++j) */
  
  /* freeing archives */
  free(path);
  for (i = 0; i < TOT_ARCH; i++)
    free_archive(arch[i],N[i]);
}

/* print parameters */
void pr_parameters(void)
{
  FILE *f;
  char dir[MAXSTR] = {'\0'};
  const char *path = 0;
  int i = 0;
  
  
  if (get_parameter("Diagnostics"))
    path = PgetsEZ("Diagnostics");
  else
    path = Pgets("top_directory");
  
  sprintf(dir,"%s/parameters.out",path);
  f = Fopen(dir,"a");
  
  fprintf(f,SECTION"Parameters"SECTION"\n");
  
  while (parameters_global[i] != 0)
  {
    fprintf(f,"%s = %s\n",parameters_global[i]->lv,parameters_global[i]->rv);
    i++;
  }
  
  Fclose(f);
}

/* print coords */
void pr_coords(const Grid_T *const grid)
{
  FILE *f;
  char dir[MAXSTR]={'\0'}, *path;
  const char *path_par;
  Uint i = 0;
  
  if (get_parameter("Diagnostics"))
    path_par = PgetsEZ("Diagnostics");
  else
    path_par = PgetsEZ("top_directory");
    
  path = make_directory(path_par,"Patches");
  
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    Uint U = countf(patch->node);
    Uint l;
    
    sprintf(dir,"%s/%s.patch",path,patch->name);
    f = Fopen(dir,"w");
    
    for (l = 0; l < U; l++)
      fprintf(f,"%f %f %f\n",
        patch->node[l]->x[0],
          patch->node[l]->x[1],
            patch->node[l]->x[2]);
    
    Fclose(f);
  }
  free(path);
}

/* print the difference of two given name of fields at each node 
// on the whole grid.
*/
void pr_field_difference(const Grid_T *const grid,const char *const fld1,const char *const fld2)
{
  FILE *file1,*file2;
  char dir[MAXSTR]={'\0'}, *path;
  const char *path_par;
  Uint l,i,R;
  
  path_par = Pgets("top_directory");
  path = make_directory(path_par,"Fields");
  
  sprintf(dir,"%s/%s-%s.grid",path,fld2,fld1);
  file1 = Fopen(dir,"w");
  
  fprintf(file1,"# node %s %s |%s-%s|\n",fld2,fld1,fld2,fld1);
  R = 0;
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    const double *f1 = patch->fields[Ind(fld1)]->v;
    const double *f2 = patch->fields[Ind(fld2)]->v;
    Uint U = patch->nn;
    
    sprintf(dir,"%s/%s-%s.%s",path,fld2,fld1,patch->name);
    file2 = Fopen(dir,"w");
    fprintf(file2,"# node %s %s |%s-%s|\n",fld2,fld1,fld2,fld1);
    
    for (l = 0; l < U; l++)
    {
      fprintf(file1,"%u %0.15f %0.15f %0.15f\n",l+R,f2[l],f1[l],fabs(f2[l]-f1[l]));
      fprintf(file2,"%u %0.15f %0.15f %0.15f\n",l,f2[l],f1[l],fabs(f2[l]-f1[l]));
    }
      
    R += U;
    Fclose(file2);
  }
  
  Fclose(file1);
  free(path);
}


/* adding s1 and s2 to archive to make their name 
// for printing name purposes. note: if one if s1 or s2 is null,
// the none null one will be written in A[i].s1.
*/
static void add_to_archive(struct Archive_S **const arch,SubFace_T *const s1,SubFace_T *const s2,Uint *const n,const char *const desc)
{
  struct Archive_S *A = (*arch);
  char str[MAXSTR] = {'\0'};
  const Uint N = *n;
  
  if (s1 == 0 && s2 == 0) 
    return;
  
  else if (s1 && s2)
  {
    Uint j;
    
    /* skiping for repetitive structures */
    for (j = 0; j < N; j++)
    {
      if ( (A[j].s1 == s1 && A[j].s2 == s2) ||
           (A[j].s1 == s2 && A[j].s2 == s1)   ) 
        return;
    }
    
    A = realloc(A,(N+1)*sizeof(*A));
    IsNull(A);
    
    A[N].s1 = s1;
    A[N].s2 = s2;
  
    str[0] = '\0';
    sprintf(str,"%s:%u:a(p:%s,f:%u,sf:%u)",desc,N,s1->patch->name,s1->face,s1->sn);
    A[N].n1 = dup_s(str);
    
    str[0] = '\0';
    sprintf(str,"%s:%u:b(p:%s,f:%u,sf:%u)",desc,N,s2->patch->name,s2->face,s2->sn);
    A[N].n2 = dup_s(str);
  }
  else if (s1 && !s2)
  {
    A = realloc(A,(N+1)*sizeof(*A));
    IsNull(A);
    
    A[N].s1 = s1;
    A[N].s2 = 0;
  
    str[0] = '\0';
    sprintf(str,"%s:%u:(p:%s,f:%u,sf:%u)",desc,N,s1->patch->name,s1->face,s1->sn);
    A[N].n1 = dup_s(str);
  }
  else if (!s1 && s2)
  {
    A = realloc(A,(N+1)*sizeof(*A));
    IsNull(A);
    
    A[N].s1 = s2;
    A[N].s2 = 0;
  
    str[0] = '\0';
    sprintf(str,"%s:%u(p:%s,f:%u,sf:%u)",desc,N,s2->patch->name,s2->face,s2->sn);
    A[N].n1 = dup_s(str);
  }
  
  (*arch) = A;
  ++(*n);
}

/* freeing archive */
static void free_archive(struct Archive_S *arch,const Uint N)
{
  Uint i;
  
  for(i = 0; i < N; ++i)
  {
    free(arch[i].n1);
    
    if (arch[i].s2)
      free(arch[i].n2);
  }
  
  free(arch);
}

/* print derivatives numc[#]-anac[#] versus node # for each given patch.
// ->return valuenn: the absolute value of the greater difference of  numc[#]-anac[#]. */
double pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const Patch_T *const patch,const char *const prefix)
{
  FILE *f;
  double max = 0;/* greater difference */
  char file_name[MAXSTR];
  Uint nn;
  Uint p;
  
  if (!numc)
    Error0("There is no numeric value.\n");
  
  if (!anac)
    Error0("There is no analytic value.\n");
  
  nn = total_nodes_patch(patch);
  sprintf(file_name,"%s.%s",prefix,patch->name);
  f = Fopen(file_name,"w");
  
  fprintf(f,"#Node (df/d?|N-df/d?|A) df/d?|N df/d?|A i j k x y z:\n");
  for (p = 0; p < nn; ++p)
  {
    Uint i1,j1,k1;
    double diff = numc[p]-anac[p];
    double abs_diff = fabs(diff);
    ijk_to_i_j_k(p,patch->n,&i1,&j1,&k1);
    fprintf(f,"%u %g %g %g %u %u %u %g %g %g\n",p,diff,numc[p],anac[p],i1,j1,k1,x_(p),y_(p),z_(p));
    
    if (GRT(abs_diff,max))
      max = abs_diff;
  }
  Fclose(f);
  
  return max;
}

/* printing the given matrix with specified name in "Matrix" folder.
*/
void pr_matrix(const Matrix_T *const M,const char *const name)
{
  char *folder = open_folder("Matrix");
  char path[MAXSTR];
  FILE *file;
  long r,c;
  
  sprintf(path,"%s/%s",folder,name);
  file = Fopen(path,"w");
  
  fprintf(file,"#row = %ld ,#column = %ld\n",M->row,M->col);
  
  if (!M->row || !M->col)
  {
    Fclose(file);
  }
  else if (M->reg_f)
  {
    double **const m = M->reg->A;
    
    for (r = 0; r < M->row; ++r)
    {
      for (c = 0; c < M->col; ++c)
        fprintf(file,"%g ",m[r][c]);
      fprintf(file,"\n");
    }
  }
  else if (M->tri_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (M->ccs_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (M->crs_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (M->tri_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (M->ccs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else if (M->crs_l_f)
  {
    Error0(INCOMPLETE_FUNC);
  }
  else
    Error0("No matrix format is defined for this given matrix.\n");

  Fclose(file);
  free(folder);
}
