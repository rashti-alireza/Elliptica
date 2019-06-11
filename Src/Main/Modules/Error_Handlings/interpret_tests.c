/*
// Alireza Rashti
// July 2018
*/

#include "interpret_tests.h"

/* check if the status for result of test. */
void check_test_result(const int status)
{
  switch(status)
  {
    case TEST_SUCCESSFUL:
      printf("PASSED :)\n");
      break;
    case TEST_UNSUCCESSFUL:
      printf("FAILED :(\n");
      break;
  }
}

