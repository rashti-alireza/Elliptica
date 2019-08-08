/* defining macro such that only gets the field value if the patch covers the BH horizon. */
#define GET_FIELD_IF_ON_HORIZON(x) \
 double *xNAME = 0;/* so it gets segfault if mistakenly the value is acquired. */\
 if (regex_search("grid[[:digit:]]+_right_BH_surrounding_.+",patch->name))\
   xNAME = patch->pool[Ind(#xNAME)]->v;


void *eq_psi(void *vp1,void *vp2);
void *bc_psi(void *vp1,void *vp2);
void *jacobian_eq_psi(void *vp1,void *vp2);
void *jacobian_bc_psi(void *vp1,void *vp2);

/*
void *eq_eta(void *vp1,void *vp2);
void *bc_eta(void *vp1,void *vp2);
void *jacobian_eq_eta(void *vp1,void *vp2);
void *jacobian_bc_eta(void *vp1,void *vp2);

void *eq_Beta_U0(void *vp1,void *vp2);
void *bc_Beta_U0(void *vp1,void *vp2);
void *jacobian_eq_Beta_U0(void *vp1,void *vp2);
void *jacobian_bc_Beta_U0(void *vp1,void *vp2);

void *eq_Beta_U1(void *vp1,void *vp2);
void *bc_Beta_U1(void *vp1,void *vp2);
void *jacobian_eq_Beta_U1(void *vp1,void *vp2);
void *jacobian_bc_Beta_U1(void *vp1,void *vp2);

void *eq_Beta_U2(void *vp1,void *vp2);
void *bc_Beta_U2(void *vp1,void *vp2);
void *jacobian_eq_Beta_U2(void *vp1,void *vp2);
void *jacobian_bc_Beta_U2(void *vp1,void *vp2);

void *eq_phi(void *vp1,void *vp2);
void *bc_phi(void *vp1,void *vp2);
void *jacobian_eq_phi(void *vp1,void *vp2);
void *jacobian_bc_phi(void *vp1,void *vp2);


*/