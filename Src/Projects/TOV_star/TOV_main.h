#include "core_lib.h"
#include "physics_lib.h"
#include "TOV_lib.h"

int TOV_star(void *vp);
TOV_T *TOV_solution(TOV_T *const TOV);
TOV_T *TOV_init(void);
void TOV_free(TOV_T *TOV);
