/*
// Alireza Rashti
// August 2019
*/

#include "bbn_study_bbn.h"

/* study and analyze initial data. whenever it is called, 
// it increments the solving_iteration_number by 1. */
void bbn_study_initial_data(Grid_T *const grid)
{
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
  
  Pseti("solving_iteration_number",++solving_iter);
  
  /* Observables */
  Observable_T *obs = init_observable(grid);
  obs->quantity = "ADM_momentums";
  plan_observable(obs);
  double P_ADM[3];
  double J_ADM[3];
  
  P_ADM[0] = obs->Px_ADM(obs);
  P_ADM[1] = obs->Py_ADM(obs);
  P_ADM[2] = obs->Pz_ADM(obs);
  
  J_ADM[0] = obs->Jx_ADM(obs);
  J_ADM[1] = obs->Jy_ADM(obs);
  J_ADM[2] = obs->Jz_ADM(obs);
  
  printf("ADM momentums:\n");
  printf("P_ADM = (%e,%e,%e).\n",P_ADM[0],P_ADM[1],P_ADM[2]);
  printf("J_ADM = (%e,%e,%e).\n",J_ADM[0],J_ADM[1],J_ADM[2]);
  
  free_observable(obs);

  printf("} Studying initial data for binary BH and NS ==> Done.\n");
  pr_clock();
  pr_line_custom('=');
}

/* print residual norms L2, L1 and L_inf of the specified fields. */
void bbn_print_residual_norms(Grid_T *const grid,const unsigned iteration, const char *const folder)
{
  /* list of the fields to be printed out */
  const char *f[] = {"ham_constraint",
                     "mom_constraint_U0",
                     "mom_constraint_U1",
                     "mom_constraint_U2",
                     "ham_constraint_2nd",
                     "mom_constraint_2nd_U0",
                     "mom_constraint_2nd_U1",
                     "mom_constraint_2nd_U2",
                     "B0_U0_residual",
                     "B0_U1_residual",
                     "B0_U2_residual",
                     "psi_residual",
                     "eta_residual",
                     "phi_residual",
                     0};
  unsigned i,p;
                     
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
      
      fclose(file_Linf);
      fclose(file_L1);
      fclose(file_L2);
    }
  }
  
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
