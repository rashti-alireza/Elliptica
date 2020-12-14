#include "adm_header.h"
#include "physics_freedata_lib.h"

/* handy for comparison */
#define Is_Different(x,pr) \
  if (!EQL((x),0.))\
  {\
    if (pr)\
      printf("%s = %0.1e\n",#x,(x));\
    max = (max < x ? x : max);\
  }

void adm_doctest_AConfIJ(Physics_T *const phys);
static void KerrSchild_analytic_adm_Kij_comparison(Physics_T *const phys);
