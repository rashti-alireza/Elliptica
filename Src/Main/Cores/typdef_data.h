/*
// Alireza Rashti
// May 2018
*/

/* common flags*/
typedef enum FLAG_T
{
  NUMERIC,
  LITERAL,
  YES
}Flag_T;

/* print flags */
typedef enum PRINT_T
{
  PRINT_PARAMETER
}Print_T;
/* parameters */
typedef struct PARAMETER_T
{
  /* syntax is expected to be lv = rv */
  char *lv; //letf value
  char *rv;// right value
}Parameter_T;

