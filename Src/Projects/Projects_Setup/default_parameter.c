/*
// Alireza Rashti
// June 2018
*/

#include "default_parameter.h"

/* setting the default value of parameters which are needed 
// if the user needs to make sure that the parameter must have default
// value, they need to add the parameter in this function, 
// like the example.
*/
void set_default_parameter(void)
{
  /* set the right hand side to the left hand side 
  // as the default value
  */
  set_default("output_directory_name",inputfile_name_global); 
  set_default("Fourier_Transformation_Method","RFT"); 
  
  
}

static void set_default(const char *const lhs,const char *const rhs)
{
  const char *v;
  Parameter_T *par;
  
  par = get_parameter(lhs);
  if (par == 0)
    add_parameter(lhs,rhs);
  else
  {
    v = PgetsEZ(lhs);
    if (v == 0)
      par->rv = dup_s(rhs);
    else if (v[0] == '\0' || strcmp_i(v,"default"))
    {
      free(par->rv);
      par->rv = dup_s(rhs);
    }
  }
}
