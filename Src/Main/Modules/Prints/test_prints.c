/*
// Alireza Rashti
// June 2018
*/

#include "test_prints.h"
#define MAXSTR 400
/* printing different quantities for test */
/* ->return value: if print option is on return 1 otherwise 0 */
int test_print(const Print_T f)
{
  char *on; 
  
  switch(f)
  {
    case PRINT_PARAMETERS:
      on = get_parameter_value_S("print_parameters",0);
      if (on == 0) return 0;
      if (strcmp_i(on,"yes")|| strcmp_i(on,"y"))
        return 1;
      break;
    case PRINT_COORDS:
      on = get_parameter_value_S("print_coords",0);
      if (on == 0) return 0;
      if (strcmp_i(on,"yes")|| strcmp_i(on,"y"))
        return 1;
      break;
    case PRINT_INTERFACES:
      on = get_parameter_value_S("print_interfaces",0);
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
  char str[MAXSTR]={'\0'}, *path;
  Interface_T **face;
  Node_T **node;
  SubFace_T *subf,*subf2;
  Flag_T flg;
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
  unsigned N[TOT_ARCH];
  unsigned pa,fc,sf,i,j,n,save;
  
  /* initializing */
  for (i = 0; i < TOT_ARCH; ++i)
  {
    arch[i] = 0;
    N[i] = 0;
  }
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  
  path = dup_s(path);
  make_directory(&path,"InterfaceInfo",YES);
  
  str[0] = '\0';
  sprintf(str,"%s/interface_pairings.rm",path);
  f = fopen(str,"w");
  pointerEr(f);
  
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
  
  fclose(f); 
  
  /* printing interfacess */
  for (n = 0; n < TOT_ARCH; ++n)
  {
    for (i = 0; i < N[n]; ++i)
    {
      if (arch[n][i].s1 && arch[n][i].s2)
      {
        str[0] = '\0';
        sprintf(str,"%s/%s.paired",path,arch[n][i].n1);
        f = fopen(str,"w");
        pointerEr(f);
        fprintf(f,"#%s\n",arch[n][i].s1->flags_str);
        
        node = arch[n][i].s1->patch->node;
        
        for (j = 0; j < arch[n][i].s1->np; ++j)
        {
          double *x = node[arch[n][i].s1->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        fclose(f);
        
        str[0] = '\0';
        sprintf(str,"%s/%s.paired",path,arch[n][i].n2);
        f = fopen(str,"w");
        pointerEr(f);
        fprintf(f,"#%s\n",arch[n][i].s2->flags_str);
        
        node = arch[n][i].s2->patch->node;
        
        for (j = 0; j < arch[n][i].s2->np; ++j)
        {
          double *x = node[arch[n][i].s2->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        fclose(f);
      }/* if (arch[n][i].s1 && arch[n][i].s2) */
      else
      {
        str[0] = '\0';
        sprintf(str,"%s/%s.single",path,arch[n][i].n1);
        f = fopen(str,"w");
        pointerEr(f);
        fprintf(f,"#%s\n",arch[n][i].s1->flags_str);
        
        node = arch[n][i].s1->patch->node;
        
        for (j = 0; j < arch[n][i].s1->np; ++j)
        {
          double *x = node[arch[n][i].s1->id[j]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        fclose(f);
      
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
  char dir[MAXSTR]={'\0'}, *path;
  int i = 0;
  Flag_T flg;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  sprintf(dir,"%s/parameters.out",path);
  f = fopen(dir,"w");
  pointerEr(f);
  
  fprintf(f,SECTION"Parameters"SECTION"\n");
  
  while (parameters_global[i] != 0)
  {
    fprintf(f,"%s = %s\n",parameters_global[i]->lv,parameters_global[i]->rv);
    i++;
  }
  
  fclose(f);
}

/* print coords */
void pr_coords(const Grid_T *const grid)
{
  FILE *f;
  char dir[MAXSTR]={'\0'}, *path;
  unsigned i = 0;
  Flag_T flg;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    unsigned U = countf(patch->node);
    unsigned l;
    
    sprintf(dir,"%s/%s.patch",path,patch->name);
    f = fopen(dir,"w");
    pointerEr(f);
    
    for (l = 0; l < U; l++)
      fprintf(f,"%f %f %f\n",
        patch->node[l]->x[0],
          patch->node[l]->x[1],
            patch->node[l]->x[2]);
    
    fclose(f);
  }
  
}

/* adding s1 and s2 to archive to make their name 
// for printing name purposes. note: if one if s1 or s2 is null,
// the none null one will be written in A[i].s1.
*/
static void add_to_archive(struct Archive_S **const arch,SubFace_T *const s1,SubFace_T *const s2,unsigned *const n,const char *const desc)
{
  struct Archive_S *A = (*arch);
  char str[MAXSTR] = {'\0'};
  const unsigned N = *n;
  
  if (s1 == 0 && s2 == 0) 
    return;
  
  else if (s1 && s2)
  {
    unsigned j;
    
    /* skiping for repetitive structures */
    for (j = 0; j < N; j++)
    {
      if ( (A[j].s1 == s1 && A[j].s2 == s2) ||
           (A[j].s1 == s2 && A[j].s2 == s1)   ) 
        return;
    }
    
    A = realloc(A,(N+1)*sizeof(*A));
    pointerEr(A);
    
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
    pointerEr(A);
    
    A[N].s1 = s1;
    A[N].s2 = 0;
  
    str[0] = '\0';
    sprintf(str,"%s:%u:(p:%s,f:%u,sf:%u)",desc,N,s1->patch->name,s1->face,s1->sn);
    A[N].n1 = dup_s(str);
  }
  else if (!s1 && s2)
  {
    A = realloc(A,(N+1)*sizeof(*A));
    pointerEr(A);
    
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
static void free_archive(struct Archive_S *arch,const unsigned N)
{
  unsigned i;
  
  for(i = 0; i < N; ++i)
  {
    free(arch[i].n1);
    
    if (arch[i].s2)
      free(arch[i].n2);
  }
  
  free(arch);
}

/* print derivatives numc[#]-anac[#] versus node # for each given patch.
// ->return valuenn: if the mentioned difference greater than 
// the given tolerance NO, YES otherwise.
*/
Flag_T pr_derivatives_DiffByNode(const double *const numc, const double *const anac,const Patch_T *const patch,const char *const prefix)
{
  FILE *f;
  char file_name[MAXSTR];
  unsigned nn;
  unsigned p;
  double tol = ROUND_OFF_ERR;
  const char *par = get_parameter_value_S("test_derivative",0);
  char *save,*tol_s = dup_s(par);
  Flag_T flg = YES;
  
  if (!numc)
    abortEr("There is no numeric value.\n");
  
  if (!anac)
    abortEr("There is no analytic value.\n");
  
  if (strchr(tol_s,COMMA))
  {
    par = tok_s(tol_s,COMMA,&save);
    tol = atof(save);
  }
  
  nn = total_nodes_patch(patch);
  sprintf(file_name,"%s.%s",prefix,patch->name);
  f = fopen(file_name,"w");
  pointerEr(f);
  
  fprintf(f,"#Node (df/d?|N-df/d?|A) df/d?|N df/d?|A i j k x y z:\n");
  for (p = 0; p < nn; ++p)
  {
    unsigned i1,j1,k1;
    double diff = numc[p]-anac[p];
    IJK(p,patch->n,&i1,&j1,&k1);
    fprintf(f,"%u %f %f %f %u %u %u %f %f %f\n",p,diff,numc[p],anac[p],i1,j1,k1,x_(p),y_(p),z_(p));
    
    if (GRT(ABS(diff),tol))
      flg = NO;
  }
  fclose(f);
  
  free(tol_s);
  return flg;
}

