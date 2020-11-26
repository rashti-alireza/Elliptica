#ifndef physics_StressEnergyTensor_LIB_H
#define physics_StressEnergyTensor_LIB_H

/* forward declaration */
struct PHYSICS_T;
struct GRID_T;


int Tij_main(struct PHYSICS_T *const phys);
int Tij_mount(struct GRID_T *const grid);


#endif


