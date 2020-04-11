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
// print_fields_4d = yes,Format:HDF5,{(V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2)}
// # as one can see the vector quantities determined by parenthesis 
// # and all of the desired fields need to be put in curly bracket
//
// Pr_Field_T *pr  = init_PrField(grid);
// pr->folder      = "folder_path";
// pr->par         = "print_fields_4d";
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
#define DL_OC '{'
#define DL_CC '}'
#define DL_OP '('
#define DL_CP ')'
#define DL_C ','
/* error msg */
#define FORMAT_ER_PAR "Incorrect parameter format.\n"\
    "Print_Fields_4D must be written like:\n"\
      "Print_Fields_4D = yes: format:..,{field1,field2,...,(fx,fy,fz),...}\n"

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and prints the whole grid and fields in the specified folder. */
void pr_fields(Pr_Field_T *const pr)
{
  const char *par4 = PgetsEZ(pr->par);/* print_fields_4d */
  const char *par3 = PgetsEZ(pr->par);/* print_fields_3d */
  const char *par2 = PgetsEZ(pr->par);/* print_fields_2d */
  
  /* 4d prints, i.e. field versus a,b,c */
  if (strstr_i(par4,"yes"))
  {
    read_parameter_4d(par4,pr);
    if (strstr_i(par4,"Format:HDF5"))
    {
      pr_hdf5_silo(pr);
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
/* initiating a Pr_Field_T for printing with given grid.
// ->return value: Pr_Field_T */
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
  _free(info);
  _free(pr);
}

/* it reads the parameter to find out what and how it's going to be print.
// given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on.
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {(field1)vs(a,b,c)|coord}. */
void read_parameter_4d(const char *const par,Pr_Field_T *const pr)
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

