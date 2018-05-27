/*
// Alireza Rashti
// May 2018
*/

/* common flags*/
typedef enum FLAG_T
{
  NUMERIC,
  LITERAL
}Flag_T;

/* parameters */
typedef struct PARAMETER_T
{
  /* syntax is expected to be lv = rv */
  char *lv; //letf value
  char *rv;// right value
}Parameter_T;

