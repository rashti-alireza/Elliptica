/*
// Alireza Rashti
// July 2018
*/

/* Print Fields:
// note: in Silo language, I'm using curvilinear format for both mesh and data (fields)
// regardless of the patch is Cartesian or not.
// note: data and mesh in Silo must be written in column major order,
//       otherwise for inhomogeneous resolutions you get a scrambled mesh.
//
// usage examples:
// ===============
// # parameter that is determined in input file is like:
// silo_output_3d = (V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2)
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

#include "pr_for_fields.h"

/* DeLimits */
#define DL_OP '('
#define DL_CP ')'
#define DL_C ','

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and prints the whole grid and fields in the specified folder. */
void pr_fields(Pr_Field_T *const pr)
{
  if (!pr)
    return;
  
  /* using silo */  
  if (PgetsEZ("silo_output_3d"))
  {
    parse_parameter_3d(Pgets("silo_output_3d"),pr);
    pr_hdf5_silo(pr);
  }
  
}
/* initiating a Pr_Field_T for printing with given grid.
// ->return value: Pr_Field_T */
Pr_Field_T *init_PrField(const Grid_T *const grid)
{
  Pr_Field_T *pr = calloc(1,sizeof(*pr));
  IsNull(pr);
  pr->grid = grid;
  
  return pr;
}

/* freeing Pr_Field_T */
void free_PrField(Pr_Field_T *pr)
{
  struct Info_S *info = pr->group;
  Uint i;
  
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
  Free(info);
  Free(pr);
}

/* it reads the parameter to find out what and how it's going to be print.
// given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on.
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {(field1)vs(a,b,c)|coord}. 
// NOTE: it supports regular expression for scalar, ex: ^dpsi_D.+ */
static void parse_parameter_3d(const char *const par,Pr_Field_T *const pr)
{
  struct Info_S *info_s = 0,*Pinfo;
  char *tok = dup_s(par);
  char *save = 0,*sub_tok = 0;
  Uint Ninfo = 0;
  
  save    = tok;
  sub_tok = tok_s(sub_tok,DL_C,&save);
  /* either => sub_tok = f1 and save = f2,f3,... 
  // or     => sub_tok = (fx,fy,fz) and save = ... */
  if (sub_tok == 0)
    fprintf(stderr,"No field in parameter file is specified for 3d printing.\n");
    
  while (sub_tok)
  {
    char *savess = 0,*ss = 0,*dump;
    
    info_s = realloc(info_s,(Ninfo+1)*sizeof(*info_s)); 
    IsNull(info_s);
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
      ss = tok_s(0,DL_CP,&savess); /* => ss = fz */
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

/* a quick print of given grid and all of the fields determined
// in the parameter file.*/
int print_fields_3D(const Grid_T *const grid,const int iteration,
                    const char *const dir)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder      = dir;
  pr->cycle       = iteration;
  pr_fields(pr);
  free_PrField(pr);

  FUNC_TOC
  return EXIT_SUCCESS;
}

/* ->: largest L2 norm.
// print norms L2, L1 and L_inf of the specified fields 
// in parameter "txt_output_1d" (supporting regular expression) */
double print_fields_1D(const Grid_T *const grid,const int iteration, 
                      const char *const folder)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  /* list of the fields to be printed out */
  char **f = 
     read_separated_items_in_string(PgetsEZ("txt_output_1d"),',');
  double largest_L2_error = 0; 
  Uint i,p;
  
  if(f)         
  for (i = 0; f[i]; ++i)
  {
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      Field_T *field = 0;
      Uint nn        = patch->nn;
      int field_ind  = _Ind(f[i]);
      Uint Nf        = 0;
      Uint *fInd     = 0;
      Uint Ufield_ind= 0;
      Uint fi;
      double Linf,L2,L1;
      FILE *file_Linf,*file_L1,*file_L2;
      char file_name_Linf[STR_LEN];
      char file_name_L1[STR_LEN];
      char file_name_L2[STR_LEN];
      char *stem = strstr(patch->name,"_");
      stem++;
      
      /* if couldn't find, maybe its a regex. */
      if (field_ind < 0)
      {
        fInd = find_field_index_with_regex(patch,f[i],&Nf);
        
        /* if nothing found. */
        if (!fInd)
          continue;
      }
      else
      {
        Ufield_ind = (Uint)field_ind;
        fInd = &Ufield_ind;
        Nf   = 1;
      }
      
      for (fi = 0; fi < Nf; ++fi)
      {
        field = patch->fields[fInd[fi]];
        if (!field->v)/* if field empty */
          continue;
          
        sprintf(file_name_Linf,"%s/%s_Linf_%s.txt",folder,field->name,stem);
        sprintf(file_name_L1,  "%s/%s_L1_%s.txt" ,folder,field->name,stem);
        sprintf(file_name_L2,  "%s/%s_L2_%s.txt" ,folder,field->name,stem);
        
        if (access(file_name_Linf,F_OK) != -1)/* if file exists */
        {
          file_Linf = Fopen(file_name_Linf,"a");
        }
        else
        {
          file_Linf = Fopen(file_name_Linf,"w");
          fprintf(file_Linf,"#iteration  %s\n",field->name);
        }
        
        if (access(file_name_L1,F_OK) != -1)
        {
          file_L1 = Fopen(file_name_L1,"a");
        }
        else
        {
          file_L1 = Fopen(file_name_L1,"w");
          fprintf(file_L1,"#iteration  %s\n",field->name);
        }
        
        if (access(file_name_L2,F_OK) != -1)
        {
          file_L2 = Fopen(file_name_L2,"a");
        }
        else
        {
          file_L2 = Fopen(file_name_L2,"w");
          fprintf(file_L2,"#iteration  %s\n",field->name);
        }
          
        Linf  = L_inf(nn,field->v);
        L2    = L2_norm(nn,field->v,0);
        L1    = L1_norm(nn,field->v,0);
        
        fprintf(file_Linf,"%-11u %0.15f\n",iteration,Linf);
        fprintf(file_L1,  "%-11u %0.15f\n",iteration,L1);
        fprintf(file_L2,  "%-11u %0.15f\n",iteration,L2);
        
        largest_L2_error = L2 > largest_L2_error ? L2 : largest_L2_error;
        
        fclose(file_Linf);
        fclose(file_L1);
        fclose(file_L2);
      }
      /* only if it is a regex found free it. */
      if (field_ind < 0)
        Free(fInd);
    }
  }
  free_2d(f);
  
  FUNC_TOC
  return largest_L2_error;
}
