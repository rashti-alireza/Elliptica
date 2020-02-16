/*
// Alireza Rashti
// June 2018
*/

#include "pr.h"

/* printing a line with a comment inside */
void pr_comment(const char *const comment)
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
  fflush(stdout);
  
}

/* print a line with a costum character*/
void pr_line_custom(const char c)
{
  int i;
  
  for (i = 0; i < MAX_LENGTH; i++)
    printf("%c",c);
  printf("\n");
  fflush(stdout);
  
}

/* print a line with a costum character*/
void pr_half_line_custom(const char c)
{
  int i;
  
  for (i = 0; i < MAX_LENGTH/2; i++)
    printf("%c",c);
  printf("\n");
  fflush(stdout);
  
}

/* printing the sepnt time form the start of abc */
void pr_clock(void)
{
  time_t now = time(0);
  double t = difftime(now,initial_time_global);
  int d,h,m,s;
  d = (int)(t/(24*60*60));
  h = (int)(t-d*(24*60*60))/(60*60);
  m = (int)(t-d*(24*60*60)-(h*60*60))/60;
  s = (int)(t-d*(24*60*60)-(h*60*60)-m*60);
  printf("\nWALL-TIME: %02dd:%02dh:%02dm:%02ds = %.0fs\n",d,h,m,s,t);
  fflush(stdout);
}

/* giving the current time in second.
// ->return current time in second.
*/
double get_time_sec(void)
{
  time_t now = time(0);
  return difftime(now,initial_time_global);
}

/* printing the amount of time spent for an event in second.
// given start time in SECOND and an appropriate message, it prints
// the difference between start and current time with the message.
*/
void pr_spent_time(const double start,const char *const event)
{
  const double end = get_time_sec();
  
  printf("Spent time for %s is \"%f\"\n",event,end-start);
}

