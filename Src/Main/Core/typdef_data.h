/*
// Alireza Rashti
// May 2018
*/

/* *******************************************
// enum:
// =====
//
// NOTE: DON'T CHANGE THE INITIALIZATION NUMBER 
// OF ENUMS. THESE CONVENTIONS ARE USED THROUGHOUT
// THE CODE. CHANGING THEM MIGHT LEAD TO SEGMENTATION
// FAULT. IF YOU WANNA ADD A NEW ENUM ADD IT BEFORE
// 'NOT_INITIALIZE' ENUM.
// *******************************************
*/

/* common flags*/
typedef enum FLAG_T
{
  UNDEFINED = 0,
  NONE,
  NO,
  YES,
  READY,
  NOT_READY,
  INUSE,
  FOUND,
  CLEAN,
  BRUTE_FORCE,
  FATAL,
  INITIALIZE,
  UP    = 0,
  DOWN  = 1,
  LEFT  = 2,
  RIGHT = 3,
  BACK  = 4,
  FRONT = 5,
  NS_T_CS,/* NS type in cubed spherical. note: NS could be any compact object like BH etc. */
  SR_T_CS,/* surrounding type in cubed spherical */
  OT_T1_CS,/* outermost type1 in cubed spherical */
  OT_T2_CS,/* outermost type2 in cubed spherical */
  NOT_INITIALIZE
}Flag_T;


/* *******************************************
// parameter:
// *******************************************
*/

/* parameters */
typedef struct PARAMETER_T
{
  /* syntax is expected to be lv = rv */
  char *lv;/* left value its name*/
  char *rv;/* right value string */
  char *rv_ip;/* right value for iterative parameter */
  double rv_double;/* right value double */
  double *rv_array;/* right value array */
  unsigned rv_n;/* right value unsigned or size of the rv_array */
  unsigned iterative;/* if this is an iterative par 1, otherwise 0. */
  unsigned double_flg;/* if this is double 1, otherwise 0 */
}Parameter_T;

/* *******************************************
// project:
// *******************************************
*/

/* functions */
typedef int ProjFunc(void);

/* projects */
typedef struct PROJECT_T
{
  char *name;/* name of project */
  char *des;/* description */
  ProjFunc *func;/* project function */
}Project_T;


