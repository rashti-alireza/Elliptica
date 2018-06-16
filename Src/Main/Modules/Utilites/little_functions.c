/*
// Alireza Rashti
// June 2018
*/

#include "little_functions.h"

/* helps you to find where litte tests start*/
void test_start(char *file,int line)
{
  printf("Test starts at\n"
    "File:%s\nLine:%d\n",file,line);
}

/* printing a line with a comment inside */
void pr_comment(char *comment)
{
  int i,d;
  
  pointerEr(comment);
  
  for (i = 0; comment[i] != '\0'; i++);
  i++;
  
  if (i%2 != 0)  i++;
  
  d = MAX_LENGTH/2-i/2;
  
  /* printing like:
  // L_SYM... comment R_SYM...
  */
  for (i = 0; i < d; i++)  printf("%c",L_SYM);
  printf(" %s ", comment);
  for (i = 0; i < d; i++)  printf("%c",R_SYM);
  printf("\n");
  
}

/* print a line */
void pr_line(void)
{
  int i;
  
  for (i = 0; i < MAX_LENGTH; i++)
    printf("%c",LINE);
  printf("\n");
  
}
