/*
// Alireza Rashti
// July 2018
*/

#include "pr_for_fields.h"
#define DELIMIT '}'
#define DELIMIT2 '{'
#define COMMA ','

/* printing all the field designated in the input file */
void pr_fields(const Grid_T *const grid)
{
  const char *par4 = get_parameter_value_S("print_fields_4d",0);
  const char *par3 = get_parameter_value_S("print_fields_3d",0);
  const char *par2 = get_parameter_value_S("print_fields_2d",0);
  char **field_name = 0, **coord[3] = {0};
  const char *path_par;
  char *folder = 0;
  unsigned nf = 0;/* number of fields */
  Flag_T flg;
  unsigned i;
  
  /* 4d prints, i.e. field versus a,b,c */
  if (strstr_i(par4,"yes"))
  {
    path_par = get_parameter_value_S("output_directory_path",&flg);
    parameterEr(flg);
    folder = dup_s(path_par);
    make_directory(&folder,"Pr_Fields_3D",YES);

    get_field_vs_coords_3d(par4,&field_name,&coord[0],&coord[1],&coord[2],&nf);
    UNUSED(grid);
    //pr_fields_on_grid(grid,folder);
    
    /* freeing */
    
    for (i = 0; i < nf; ++i)
    {
      free(field_name[i]);
      free(coord[0][i]);
      free(coord[1][i]);
      free(coord[2][i]);
    }
    free(field_name);
    free(coord[0]);
    free(coord[1]);
    free(coord[2]);
    free(folder);
  }
  /* 3d prints, i.e. ?*/
  if (strstr_i(par3,"yes"))
  {
    abortEr(INCOMPLETE_FUNC);
  }
  /* 2d prints, i.e. ?*/
  if (strstr_i(par2,"yes"))
  {
    abortEr(INCOMPLETE_FUNC);
  }
  
}

/* given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on them. 
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {field_name,a,b,c}.
*/
static void get_field_vs_coords_3d(const char *const par,char ***const f_name,char ***const c0,char ***const c1,char ***const c2,unsigned *const nf)
{
  char *tok = dup_s(par);
  char **field_name, **coord0,**coord1,**coord2;
  char *save = 0,*sub_tok = 0;
  char *info;
  
  /* initializing */
  *nf = 0;
  field_name = 0;
  coord0 = 0;
  coord1 = 0;
  coord2 = 0;
  
  /* check if parameter has been written correctly */
  if (!strchr(tok,DELIMIT) || !strchr(tok,DELIMIT2))
    abortEr("Incorrect parameter format.\n"
    "Print_Fields_4D must be written like:\n"
      "Print_Fields_4D = yes: {phi,x,y,z},{alpha,a,b,c}\n");
  
  /* sub_tokx = ...{...\0 */
  sub_tok = tok_s(tok,DELIMIT,&save);
  while (sub_tok)
  {
    char *save2 = 0,*sub_tok2 = 0;
    
    field_name = realloc(field_name,(*nf+1)*sizeof(*field_name));
    pointerEr(field_name);
    coord0 = realloc(coord0,(*nf+1)*sizeof(*coord0));
    pointerEr(coord0);
    coord1 = realloc(coord1,(*nf+1)*sizeof(*coord1));
    pointerEr(coord1);
    coord2 = realloc(coord2,(*nf+1)*sizeof(*coord2));
    pointerEr(coord2);
    
    sub_tok2 = tok_s(sub_tok,DELIMIT2,&save2); /* sub_tok2 = ...,\0 */
    info = dup_s(save2);/* info = field_name,a,b,c */
    
    sub_tok2 = tok_s(info,COMMA,&save2);/* sub_tok2 = field_name */
    field_name[*nf] = dup_s(sub_tok2);
    sub_tok2 = tok_s(0,COMMA,&save2); /* sub_tok2 = a */
    coord0[*nf] = dup_s(sub_tok2);
    sub_tok2 = tok_s(0,COMMA,&save2); /* sub_tok2 = b */
    coord1[*nf] = dup_s(sub_tok2);
    sub_tok2 = tok_s(0,COMMA,&save2); /* sub_tok2 = c */
    coord2[*nf] = dup_s(sub_tok2);
    
    sub_tok = tok_s(0,DELIMIT,&save);
    free(info);
    (*nf)++;
  }
  
  *f_name = field_name;
  *c0 = coord0;
  *c1 = coord1;
  *c2 = coord2;
  
  free(tok);
}

/* printing a scalar field on the patch 
// which the given field has been defined. 
// ->return value -> EXIT_SUCCESS.*/

/*int pr_scalar_field(const Field_T *const f)
{
  
  return EXIT_SUCCESS;
}
*/