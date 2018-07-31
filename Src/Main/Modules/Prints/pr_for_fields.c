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

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and print the result in the specified folder.
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
    
    /* freeing */
    free_info_s(pr);
    
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

static void free_info_s(Pr_Field_T *const pr)
{
  struct Info_S *info = pr->vptr;
  unsigned i;
  
  for (i = 0; i < pr->nobj; ++i)
  {
    free(info[i].field);
  }
  free(info);

}

/* printing fields with HDF5 format using silo library
static void pr_fields_on_grid_HDF5_4d(Pr_Field_T *const pr)
{
  DBfile *dbfile = 0;
  unsigned pa;
  char file_name[MAX_STR_LEN];
  
  FOR_ALL_PATCHES(pa,grid)
  {
    unsigned i;
    for (i = 0; i < nf; ++i)
    {
      
    }
  }
}*/

/* it reads the parameter to find out what and how it's going to be print.
// given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on.
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {(field1,field2...)vs(a,b,c)}.
*/
static void read_parameter_4d(const char *const par,Pr_Field_T *const pr)
{
  struct Info_S *info_s = 0,*Pinfo;
  char *tok = dup_s(par);
  char *save = 0,*sub_tok = 0;
  unsigned Ninfo = 0;
  /* check if parameter has been written correctly */
  if (!strchr(tok,DL_OC) || !strchr(tok,DL_CC))
    abortEr("Incorrect parameter format.\n"
    "Print_Fields_4D must be written like:\n"
      "Print_Fields_4D = yes: format:..,{(field1,field2,...)vs(x,y,z)}...\n");
  
  /* sub_tok = (field1,field2,...)vs(x,y,z) */
  sub_tok = sub_s(tok,DL_OC,DL_CC,&save);
  while (sub_tok)
  {
    char *savef = 0,*sub_tokf = 0;
    char *savec = 0,*sub_tokc = 0;
    char *savess = 0,*ss = 0;
    
    info_s = realloc(info_s,(Ninfo+1)*sizeof(*info_s)); 
    pointerEr(info_s);
    Pinfo = &info_s[Ninfo];
    Pinfo->field = 0;
    Pinfo->nf = 0;
    
    /* getting field names */
    sub_tokf = sub_s(sub_tok,DL_OP,DL_CP,&savef);/* sub_tokf = field1,field2,... */
    ss = tok_s(sub_tokf,DL_C,&savess);/* ss = field1*/
    while (ss)
    {
      Pinfo->field = 
        realloc(Pinfo->field,(Pinfo->nf+1)*sizeof(*Pinfo->field));
      pointerEr(Pinfo->field);
      
      sprintf(Pinfo->field[Pinfo->nf],ss);
      
      ss = tok_s(0,DL_C,&savess);
      Pinfo->nf++;
    }
    
    /* getting axis */
    sub_tokc = sub_s(savef,DL_OP,DL_CP,&savec);/* sub_tokc = a,b,c */
    ss = tok_s(sub_tokc,DL_C,&savess);/* ss = a*/
    sprintf(Pinfo->axis[0],ss);
    ss = tok_s(0,DL_C,&savess); /*ss = b*/
    sprintf(Pinfo->axis[1],ss);
    ss = tok_s(0,DL_C,&savess); /*ss = c*/
    sprintf(Pinfo->axis[2],ss);
    
    sub_tok = sub_s(0,DL_OC,DL_CC,&save);
    Ninfo++;
  }
  
  pr->vptr = info_s;
  pr->nobj = Ninfo;
  
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
  free(pr);
}
