#ifndef physics_StressEnergyTensor_LIB_H
#define physics_StressEnergyTensor_LIB_H

/* forward declaration */
struct PHYSICS_T;

int Tij_tune(struct PHYSICS_T *const phys);
int Tij_mount(Grid_T *const grid);


#endif


