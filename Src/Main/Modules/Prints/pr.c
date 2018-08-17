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
  
}

/* printing the sepnt time form the start of abc */
void pr_clock(void)
{
  time_t now = time(0);
  
  printf("Clock in Seconds: \"%f\"\n"
          "Clock in Minuets: \"%f\"\n",
            difftime(now,initial_time_global),
            difftime(now,initial_time_global)/60);
  
}

/* giving the current time in second*/
double get_time(void)
{
  time_t now = time(0);
  return difftime(now,initial_time_global)/60;
}

/* printing the amount of time spent for an event.
// given start time and an appropriate message, it prints
// the difference between start and current time with the message.
*/
void pr_spent_time(const double start,const char *const event)
{
  const double end = get_time();
  
  printf("Spent time for %s is \"%f\"\n",event,end-start);
}

