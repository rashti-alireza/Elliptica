/* some shared objects to be used for field print functions. */

/* this is an internal structure for pr_for_files functions */
struct Info_S
{
  char *field;
  char *comp[3];/* components name */
  Uint vec_flg:1;/* 1 if it is vector, 0 otherwise */
};

