/* defining macro such that only gets the field value if the patch covers the BH horizon. */
#define GET_FIELD_IF_ON_HORIZON(xNAME) \
 double *xNAME = 0;/* so it gets segfault if mistakenly the value is acquired. */\
 if (IsItHorizonPatch(patch))\
   xNAME = patch->pool[Ind(#xNAME)]->v;

/* defining macro such that only gets the field value if the patch covers the NS. */
#define GET_FIELD_IF_IN_NS(xNAME) \
 double *xNAME = 0;/* so it gets segfault if mistakenly the value is acquired. */\
 if (IsItNSPatch(patch))\
   xNAME = patch->pool[Ind(#xNAME)]->v;


void *eq_psi(void *vp1,void *vp2);
void *bc_psi(void *vp1,void *vp2);
void *jacobian_eq_psi(void *vp1,void *vp2);
void *jacobian_bc_psi(void *vp1,void *vp2);


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

