/*
// Alireza Rashti
// June 2018
*/

#include "test_prints.h"

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
  char str[1000]={'\0'}, *path;
  Interface_T **face;
  Node_T **node;
  SubFace_T *subf,*subf2;
  struct Archive_S *arch;
  Flag_T flg;
  unsigned pa,fc,sf,i,j,Narch;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  
  str[0] = '\0';
  sprintf(str,"%s/interface_pairings.out",path);
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
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn);
          fprintf(f,"%s",str);
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
    
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
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn);
          fprintf(f,"%s",str);
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
 
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
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn);
          fprintf(f,"%s",str);
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
 
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
        
        if (subf->touch == 1 && subf->copy == 0)
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n"
          ">>>flags:%s\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
 
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
        
        if (subf->touch == 1 && subf->copy == 1)
        {
          str[0] = '\0';
          sprintf(str,"o. patch:%u (%s),face:%u,subface:%u\n"
          ">>>flags:%s\n",
            subf->patch->pn,subf->patch->name,subf->face,subf->sn,subf->flags_str);
          fprintf(f,"%s",str);
          flg = FOUND;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
  
  fflush(f);
  
 /* printing out pairing information: pairing*/
  fprintf(f,"###pairing interfaces are:\n");
  flg = NONE;
  Narch = 0;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->touch == 1)
        {
          subf2 = get_paired_subface(subf);
          
          str[0] = '\0';
          sprintf(str,"|->p%u:a{ %s }\n",Narch,subf->flags_str);
          fprintf(f,"%s",str);
          str[0] = '\0';
          sprintf(str,"|->p%u:b{ %s }\n",Narch,subf2->flags_str);
          fprintf(f,"%s",str);
          
          fprintf(f,"%s\n",PR_LINE);
          
          add_to_archive(&arch,subf,subf2,Narch);
          
          flg = FOUND;
          ++Narch;
        }
       
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  if (flg == NONE)
    fprintf(f,"Nothing.\n");
 
  fclose(f); 
  
  /* printing interfacess' points: paired */
  for (i = 0; i < Narch; ++i)
  {
    str[0] = '\0';
    sprintf(str,"%s/%s.paired",path,arch[i].n1);
    f = fopen(str,"w");
    pointerEr(f);
    fprintf(f,"#%s\n",arch[i].s1->flags_str);
    
    node = arch[i].s1->patch->node;
    
    for (j = 0; j < arch[i].s1->np; ++j)
    {
      double *x = node[arch[i].s1->id[j]]->x;
      fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
    }
    fclose(f);
    
    str[0] = '\0';
    sprintf(str,"%s/%s.paired",path,arch[i].n2);
    f = fopen(str,"w");
    pointerEr(f);
    fprintf(f,"#%s\n",arch[i].s2->flags_str);
    
    node = arch[i].s2->patch->node;
    
    for (j = 0; j < arch[i].s2->np; ++j)
    {
      double *x = node[arch[i].s2->id[j]]->x;
      fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
    }
    fclose(f);
  }
  
  /* printing out the interfaces with outer-boundary*/
  unsigned NouterB = 0;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    node = grid->patch[pa]->node;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->outerB == 1)
        {
          
          str[0] = '\0';
          sprintf(str,"%s/outerB:%u(p:%u,f:%u,sf:%u)",path,NouterB,pa,fc,sf);
          f = fopen(str,"w");
          pointerEr(f);
          fprintf(f,"#%s\n",subf->flags_str);
          
          for (i = 0; i < subf->np; ++i)
          {
            double *x = node[subf->id[i]]->x;
            fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
          }
          
          NouterB++;
          fclose(f);
        }/* end of if (subf->outerB == 1) */
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  
  /* printing out the interfaces with inner-boundary*/
  unsigned NinnerB = 0;
  FOR_ALL(pa,grid->patch)
  {
    face = grid->patch[pa]->interface;
    node = grid->patch[pa]->node;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        subf = face[fc]->subface[sf];
        
        if (subf->innerB == 1)
        {
          
          str[0] = '\0';
          sprintf(str,"%s/innerB:%u(p:%u,f:%u,sf:%u)",path,NinnerB,pa,fc,sf);
          f = fopen(str,"w");
          pointerEr(f);
          fprintf(f,"#%s\n",subf->flags_str);
          
          for (i = 0; i < subf->np; ++i)
          {
            double *x = node[subf->id[i]]->x;
            fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
          }
          
          NinnerB++;
          fclose(f);
        }/* end of if (subf->outerB == 1) */
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */
  
  /* freeing archive */
  free_archive(arch,Narch);
}

/* print parameters */
void pr_parameters(void)
{
  FILE *f;
  char dir[1000]={'\0'}, *path;
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
  char dir[1000]={'\0'}, *path;
  unsigned i = 0;
  Flag_T flg;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  FOR_ALL(i,grid->patch)
  {
    Patch_T *patch = grid->patch[i];
    unsigned U = countf(patch->node);
    unsigned l;
    
    sprintf(dir,"%s/%s.out",path,patch->name);
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
// for printing name purposes 
*/
static void add_to_archive(struct Archive_S **const arch,SubFace_T *const s1,SubFace_T *const s2,const unsigned i)
{
  struct Archive_S *A = (*arch);
  char str[200] = {'\0'};
  A = realloc(A,(i+1)*sizeof(*A));
  pointerEr(A);
  
  A[i].s1 = s1;
  A[i].s2 = s2;
  
  str[0] = '\0';
  sprintf(str,"p%u:a(p:%s,f:%u)",i,s1->patch->name,s1->face);
  A[i].n1 = dup_s(str);
  
  str[0] = '\0';
  sprintf(str,"p%u:b(p:%s,f:%u)",i,s2->patch->name,s2->face);
  A[i].n2 = dup_s(str);
  
  (*arch) = A;
}

/* freeing archive */
static void free_archive(struct Archive_S *arch,const unsigned Narch)
{
  unsigned i;
  
  for(i = 0; i < Narch; ++i)
  {
    free(arch[i].n1);
    free(arch[i].n2);
  }
  
  free(arch);
}