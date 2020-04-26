/*
// Alireza Rashti
// August 2019
*/

#include "bbn_study_bbn.h"

/* study and analyze initial data. whenever it is called, 
// it increments the solving_iteration_number by 1. */
void bbn_study_initial_data(Grid_T *const grid)
{
  if (Pgeti("STOP"))
    return;
    
  pr_clock();
  pr_line_custom('=');
  printf("{ Studying initial data for binary BH and NS ...\n");
  
  if (!grid)
  {
    printf("~> grid is empty.\n");
    printf("} Studying initial data for binary BH and NS ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
    return;
  }
  
  const char *const folder = Pgets("Diagnostics");
  int solving_iter         = PgetiEZ("solving_iteration_number");
  
  /* calculating the constraints */
  bbn_calculate_constraints_1st(grid);
  bbn_calculate_constraints_2nd(grid);
  
  /* prints */
  bbn_print_fields(grid,(unsigned)solving_iter,folder);
  bbn_print_residual_norms(grid,(unsigned)solving_iter,folder);
  bbn_print_properties(grid,(unsigned)solving_iter,folder,"a",1);
  
  Pseti("solving_iteration_number",++solving_iter);
  
  printf("} Studying initial data for binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* print the properites of the system for instance:
// mass, spin, momentum, distance and etc.
// folder        : output directory 
// open_file_mode: "a" for append and "w" for newfile
// pr_flg        : if pr_flg = 1 it also prints them at stout otherwise none. */
void bbn_print_properties(Grid_T *const grid,const unsigned iteration, const char *const folder,const char *const open_file_mode,const int pr_flg)
{
  const char *const file_name = "bbn_geometry_and_physics.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
  const char *const geometry_logo = "\n############\n"
                                    "# Geometry #\n"
                                    "############\n";
  const char *const physics_logo = "\n###########\n"
                                   "# Physics #\n"
                                   "###########\n";
  
  /* open file */
  sprintf(str,"%s/%s",folder,file_name);
  file = fopen(str,open_file_mode);
  pointerEr(file);
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %u\n",iteration);
  

  if (pr_flg)
  {
    pr_line_custom('=');
    printf("{ Geometry and physics of binary BH and NS ...\n");
  }
  
  /* { geometry */
  if (pr_flg)
    printf(geometry_logo);
    
  fprintf(file,geometry_logo);
  
  PR_PARAMETR_IN_FILE(NS_center_x)
  PR_PARAMETR_IN_FILE(NS_center_y)
  PR_PARAMETR_IN_FILE(NS_center_z)
  
  PR_PARAMETR_IN_FILE(BH_center_x)
  PR_PARAMETR_IN_FILE(BH_center_y)
  PR_PARAMETR_IN_FILE(BH_center_z)
  
  PR_PARAMETR_IN_FILE(x_CM)
  PR_PARAMETR_IN_FILE(y_CM)
  PR_PARAMETR_IN_FILE(z_CM)
  
  PR_PARAMETR_IN_FILE(r_excision)
  /* } geometry */
  
  /* { physics */
  if (pr_flg)
    printf(physics_logo);
  
  fprintf(file,physics_logo);
  
  PR_PARAMETR_IN_FILE(NS_baryonic_mass)
  PR_PARAMETR_IN_FILE(BH_irreducible_mass)
  
  PR_PARAMETR_IN_FILE(P_ADM_x)
  PR_PARAMETR_IN_FILE(P_ADM_y)
  PR_PARAMETR_IN_FILE(P_ADM_z)
  
  PR_PARAMETR_IN_FILE(J_ADM_x)
  PR_PARAMETR_IN_FILE(J_ADM_y)
  PR_PARAMETR_IN_FILE(J_ADM_z)
  PR_PARAMETR_IN_FILE(largest_L2norm_error)
  /* } physics */
  
  fprintf(file,"%s\n",LINE_STR);
  fclose(file);
  
  if (pr_flg)
  {
    printf("} Geometry and physics of binary BH and NS ==> Done.\n");
    pr_clock();
    pr_line_custom('=');
  }
  
  UNUSED(grid);
}

/* print residual norms L2, L1 and L_inf of the specified fields. */
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  /* list of the fields to be printed out */
  char **f = 
     read_separated_items_in_string(PgetsEZ("output_2d_txt"),',');
  double largest_L2_error = 0;                   
  unsigned i,p;
  
  if(f)         
  for (i = 0; f[i]; ++i)
  {
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      unsigned nn    = patch->nn;
      Field_T *field;
      int field_ind = _Ind(f[i]);
      double Linf,L2,L1;
      FILE *file_Linf,*file_L1,*file_L2;
      char file_name_Linf[1000];
      char file_name_L1[1000];
      char file_name_L2[1000];
      char *stem = strstr(patch->name,"_");
      stem++;
      
      if (field_ind < 0)
        continue;
      
      field = patch->pool[field_ind];
      if (!field->v)/* if field empty */
        continue;
      
      sprintf(file_name_Linf,"%s/%s_Linf_%s.table",folder,f[i],stem);
      sprintf(file_name_L1,  "%s/%s_L1_%s.table",folder,f[i],stem);
      sprintf(file_name_L2,  "%s/%s_L2_%s.table",folder,f[i],stem);
      
      if (access(file_name_Linf,F_OK) != -1)/* if file exists */
      {
        file_Linf = fopen(file_name_Linf,"a");
        pointerEr(file_Linf);
      }
      else
      {
        file_Linf = fopen(file_name_Linf,"w");
        pointerEr(file_Linf);
        fprintf(file_Linf,"#iteration  %s\n",f[i]);
      }
      
      if (access(file_name_L1,F_OK) != -1)
      {
        file_L1 = fopen(file_name_L1,"a");
        pointerEr(file_L1);
      }
      else
      {
        file_L1 = fopen(file_name_L1,"w");
        pointerEr(file_L1);
        fprintf(file_L1,"#iteration  %s\n",f[i]);
      }
      
      if (access(file_name_L2,F_OK) != -1)
      {
        file_L2 = fopen(file_name_L2,"a");
        pointerEr(file_L2);
      }
      else
      {
        file_L2 = fopen(file_name_L2,"w");
        pointerEr(file_L2);
        fprintf(file_L2,"#iteration  %s\n",f[i]);
      }
        
      Linf  = L_inf(nn,field->v);
      L2    = L2_norm(nn,field->v,0);
      L1    = L1_norm(nn,field->v,0);
      
      fprintf(file_Linf,"%-11u %0.15f\n",iteration,Linf);
      fprintf(file_L1,  "%-11u %0.15f\n",iteration,L1);
      fprintf(file_L2,  "%-11u %0.15f\n",iteration,L2);
      
      /* since we do all sort of approx. inside the BH don't 
      // consider this as the numerical error. */
      if (!IsItInsideBHPatch(patch))
        largest_L2_error = L2 > largest_L2_error ? L2 : largest_L2_error;
      
      fclose(file_Linf);
      fclose(file_L1);
      fclose(file_L2);
    }
  }
  free_2d(f);
  /* update numeric error */
  Psetd("largest_L2norm_error",largest_L2_error);
}

/* printing fields determined in parameter file */
void bbn_print_fields(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  pr_line_custom('=');
  printf("{ Printing Specified Fields for Binary BH and NS ...\n");

  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder = folder;
  pr->par    = "print_fields_4d";
  pr->cycle  = (int)iteration;
  pr_fields(pr);
  free_PrField(pr);
  
  printf("} Printing Specified Fields for Binary BH and NS ==> Done.\n");
  pr_line_custom('=');
}
