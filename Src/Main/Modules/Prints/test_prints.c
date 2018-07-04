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
  char dir[1000]={'\0'}, *path;
  Flag_T flg;
  unsigned pa,fc,sf,i;
  
  path = get_parameter_value_S("output_directory_path",&flg);
  parameterEr(flg);
  
  FOR_ALL(pa,grid->patch)
  {
    Interface_T **face = grid->patch[pa]->interface;
    Node_T **node = grid->patch[pa]->node;
    
    FOR_ALL(fc,face)
    {
      for (sf = 0; sf < face[fc]->ns; ++sf)
      {
        SubFace_T *subf = face[fc]->subface[sf];
        
        dir[0] = '\0';
        sprintf(dir,"%s/%s",path,subf->flags_str);
        f = fopen(dir,"w");
        pointerEr(f);
        fprintf(f,"#%s\n",subf->flags_str);
        
        for (i = 0; i < subf->np; ++i)
        {
          double *x = node[subf->id[i]]->x;
          fprintf(f,"%f %f %f\n",x[0],x[1],x[2]);
        }
        
        fclose(f);
      }
    }
  }/* end of FOR_ALL(pa,grid->patch) */

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
