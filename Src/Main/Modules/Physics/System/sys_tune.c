/*
// Alireza Rashti
// November 2020
*/


/* collection of functions to deal with tuning of the physical system */

#include "sys_tune.h"

/* tuning P_ADMs to be zero. */
int sys_tune_ADM_momenta(Physics_T *const phys)
{
  FUNC_TIC
  
  const char *const par = sysGets("P_ADM_control_method");
  char *adjust[3] = {0};
  double adm_p[3] = {0};
  double adm_j[3] = {0};
  
  printf(Pretty0"%s\n",par);
  
  parse_adjust_parameter(par,adjust);
  
  void (*P_ADM_control_0)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[0]);
  
  void (*P_ADM_control_1)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[1]);
                              
  void (*P_ADM_control_2)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[2]);
  
  
  observe(phys,"ADM(P)",Gets("Observe_ADM_P"),adm_p);
  observe(phys,"ADM(J)",Gets("Observe_ADM_J"),adm_j);
  
  printf(Pretty0"Current P_ADM = (%e,%e,%e)\n",adm_p[0],adm_p[1],adm_p[2]);
  printf(Pretty0"Current J_ADM = (%e,%e,%e)\n",adm_j[0],adm_j[1],adm_j[2]);
  
  sysSetd("Px_ADM",adm_p[0]);
  sysSetd("Py_ADM",adm_p[1]);
  sysSetd("Pz_ADM",adm_p[2]);
  
  sysSetd("Jx_ADM",adm_j[0]);
  sysSetd("Jy_ADM",adm_j[1]);
  sysSetd("Jz_ADM",adm_j[2]);
  
  if (P_ADM_control_0)
    P_ADM_control_0(phys);
  
  if (P_ADM_control_1)
    P_ADM_control_1(phys);
    
  if (P_ADM_control_2)
    P_ADM_control_2(phys);
    
  Free(adjust[0]);
  Free(adjust[1]);
  Free(adjust[2]);
  
  FUNC_TOC
  
  return EXIT_SUCCESS;
}

/* getting adjustment str, returns the relevant function. */
static fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust)
{
  fAdjustment_t *f = 0;
  
  if (!adjust)
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"none"))
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"x_CM"))
  {
    f = Py_ADM_is0_by_x_CM;
  }
  else if (strcmp_i(adjust,"y_CM"))
  {
    f = Px_ADM_is0_by_y_CM;
  }
  else
    Error0(NO_OPTION);
  
  return f;
}

/* parsing adjust parameter value and fill components consequently */
static void parse_adjust_parameter(const char *const par,char *adjust[3])
{
  if (!strstr_i(par,"adjust(") && !strstr_i(par,"none"))
    Errors("Syntax error for '%s'.\n",par);
  
  /* if it is none */  
  if (strcmp_i(par,"none"))
  {
    adjust[0] = 0;
    adjust[1] = 0;
    adjust[2] = 0;
    
    return;
  }
  
  /* parse if not none */
  char *str = dup_s(par);
  char *save;
  char *main_str = sub_s(str,'(',')',&save);
  char *tok  = tok_s(main_str,',',&save);
  
  adjust[0] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[1] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[2] = dup_s(tok);
  
  free(str);
}

/* find y_CM by demanding Px_ADM = 0 */
static void Px_ADM_is0_by_y_CM(Physics_T *const phys)
{
  double dy_CM = 0,y_CM_new;
  const double px    = sysGetd("Px_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double Omega = sysGetd("angular_velocity");
  const double y_CM0 = sysGetd("y_CM");
  const double admM  = sysGetd("ADM_mass");
  
  printf(Pretty0"adjusting Px_ADM by y_CM ...\n");
  
  /* changing center of mass */
  dy_CM    = -px/(Omega*(admM));
  y_CM_new = y_CM0+dy_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(fabs(px),dP))
  {
    sysSetd("y_CM",y_CM_new);
  }
}

/* find x_CM by demanding Py_ADM = 0 */
static void Py_ADM_is0_by_x_CM(Physics_T *const phys)
{
  double  dx_CM = 0,x_CM_new;
  const double py    = sysGetd("Py_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double Omega = sysGetd("angular_velocity");
  const double x_CM0 = sysGetd("x_CM");
  const double admM  = sysGetd("ADM_mass");
  
  printf(Pretty0"adjusting Py_ADM by x_CM ...\n");
  
  /* changing center of mass */
  dx_CM    = py/(Omega*(admM));
  x_CM_new = x_CM0+dx_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(fabs(py),dP))
  {
    sysSetd("x_CM",x_CM_new);
  }
}


