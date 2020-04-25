/*
// Alireza Rashti
// May 2018
*/

#include "parameters.h"

/* add or update a string value parameter */
void update_parameter_string(const char *const lv, const char *const rv)
{
  if (!lv || !rv)
    return;
    
  Parameter_T *par;
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    _free(par->rv);
    par->rv = dup_s(rv);
    par->double_flg = 0;
  }
  else/* since it does not exist */
    add_parameter(lv,rv);
  
}

/* add or update an integer value parameter */
void update_parameter_integer(const char *const lv, const int rv)
{
  if (!lv)
    return;
    
  Parameter_T *par;
  char str_rv[100];
  sprintf(str_rv,"%d",rv);
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    _free(par->rv);
    par->rv = dup_s(str_rv);
    par->double_flg = 0;
  }
  else/* since it does not exist */
  {
    add_parameter(lv,str_rv);
  }
  
}

/* add or update a double value parameter, 
// if print_flg != 0 it prints the changes. */
void update_parameter_double_format(const char *const lv, const double rv,const int print_flg)
{
  if (!lv)
    return;
    
  Parameter_T *par;
  char str_rv[100] = {'\0'};
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    if (print_flg)
    {
      printf("Updating Parameter:\n");
      printf("         |--> parameter     = %s\n",lv);
      printf("         |--> new value     = %+e\n",rv);
      printf("         |--> old value     = %+e\n",par->rv_double);
      printf("         |--> v_new - v_old = %+e\n",rv-par->rv_double);
    }
    sprintf(str_rv,"%15.18f",rv);
    _free(par->rv);
    
    /* NOTE:crucial to write in str format for checkpoint file purposes */
    par->rv         = dup_s(str_rv);
    par->rv_double  = rv;
    par->double_flg = 1;
  }
  else/* since it does not exist */
  {
    add_parameter_double(lv,rv);
  }
  
}

/* adding left value and right value to parameter data base 
// double format. */
void add_parameter_double(const char *const lv, const double rv)
{
  pointerEr(lv);
  
  Parameter_T *par;
  char str_rv[100] = {'\0'};
  
  par = get_parameter(lv);
  if (par)
    abortEr_s("This parameter \"%s\" has already been added!\n",lv);
    
  sprintf(str_rv,"%15.18f",rv);
  par = alloc_parameter(&parameters_global);
  par->rv = dup_s(str_rv);
  par->lv = dup_s(lv);
  par->rv_double  = rv;
  par->double_flg = 1;
  
  printf("Adding Parameter:\n");
  printf("       |--> parameter = %s\n",lv);
  printf("       |--> value     = %+g\n",rv);
}

/* adding left value and right value to parameter data base 
// array format. */
void add_parameter_array(const char *const lv, const double *const rv,const unsigned n)
{
  pointerEr(lv);
  pointerEr(rv);
  
  Parameter_T *par;
  unsigned i;
  
  par = get_parameter(lv);
  if (par)
    abortEr_s("This parameter '%s' has already been added!\n",lv);
    
  par = alloc_parameter(&parameters_global);
  par->lv = dup_s(lv);
  
  par->rv_array = alloc_double(n);
  par->rv_n = n;
  for (i = 0; i < n; ++i)
    par->rv_array[i] = rv[i];
  
}

/* update parameter array format. */
void update_parameter_array(const char *const lv, const double *const rv,const unsigned n)
{
  pointerEr(lv);
  pointerEr(rv);
  
  Parameter_T *par;
  unsigned i;
  
  par = get_parameter(lv);
  
  if (par)
  {
    if (par->iterative)
      abortEr_s("Wrong update: parameter '%s' is iterative.\n",lv);
    
    par->double_flg = 0;
      
    _free(par->rv_array);
  }
  else
  {  
    par = alloc_parameter(&parameters_global);
    par->lv = dup_s(lv);
  }
  
  par->rv_array = alloc_double(n);
  par->rv_n = n;
  for (i = 0; i < n; ++i)
    par->rv_array[i] = rv[i];
}

/* adding left value and right value to parameter data base 
// string format */
void add_parameter_string(const char *const lv, const char *const rv)
{
  add_parameter(lv,rv);
}

/* adding left value and right value to parameter data base 
// string format (the most common one) */
void add_parameter(const char *const lv, const char *const rv)
{
  pointerEr(lv);
  
  Parameter_T *par;
  
  par = get_parameter(lv);
  if (par)
    abortEr_s("This parameter \"%s\" has already been added!\n",lv);
    
  par = alloc_parameter(&parameters_global);
  par->lv = dup_s(lv);
  if (rv == 0)
    par->rv = 0;
  else if(!strcmp(rv,"") || !strcmp(rv," ") || rv[0] == '\0')
    par->rv = 0;
  else
  {
    /* study if this is an iterative parameter */
    if (strstr(rv,"->") || regex_search("\\(x[[:digit:]]+\\)",rv))
    {
      par->iterative = 1;
      /* if it has multiplicity */
      if (regex_search("\\(x[[:digit:]]+\\)",rv))
        par->rv_ip = parse_multiplicity_of_iterative_parameter(rv);
      else
        par->rv_ip = dup_s(rv);
        
      par->rv = get_n_value_str_ip(par,0);/* setting the first value of the iterative parameter */
    }
    else
    {
      par->iterative = 0;
      par->rv = dup_s(rv);
    }
  }
}

/* some iterative parameter may have multiplicity with (x?), e.g:
// par = 1->2(x2)->5(x3)->8
// this function parse it and return it in regular format, i.e:
// par = 1->2->2->5->5->5->8
// ->return value: an iterative parameter with regular format */
static char *parse_multiplicity_of_iterative_parameter(const char *const rv)
{
  if (regex_search("->$",rv))
    abortEr_s("Wrong syntax for '%s'; '->' at the end of the line.",rv);
    
  char *ret_str = 0;
  char *str = dup_s(rv);
  
  char *m_str = regex_find("\\(x[[:digit:]]+\\)",rv);
  while (m_str)
  {
    /* reference: 1->2(x2)->5(x3)->8 */
    unsigned len = (unsigned)strlen(m_str);/* (x2) => 4 */
    const char *l_str = strstr(str,m_str);/*    => (x2)->5(x3)->8 */
    const char *r_str = strlen(l_str) > len ? l_str+len: "\0";/* => ->5(x3)->8 */
    
    const char *b_str = l_str;
    while (b_str != str && *b_str != '>')
      b_str--;
    /* b_str => >2 */  
    if (*b_str == '>')
      b_str++;/* b_str => 2 */
    
    unsigned n = (unsigned)(l_str-b_str+1);  
    char *v_str = calloc(n,1);
    pointerEr(v_str);
    
    unsigned i = 0;
    while (i < n-1)
    {
      v_str[i] = b_str[i];
      i++;
    }/* =? v_str = 2 */
    v_str[n-1] = '\0';
    
    char *mult_str  = regex_find("[[:digit:]]+",m_str);
    unsigned mult   = (unsigned)atoi(mult_str);
    
    if (mult == UINT_MAX)
      abortEr_s("Wrong syntax for '%s'; negative multiplicity.\n",rv);
    if (mult == 0)
      abortEr_s("Wrong syntax for '%s'; zero multiplicity.\n",rv);
    
    unsigned n_mult = (mult-1)*2+mult*(n-1)+1;/* for each '->' 2 Byte,
                                    // for each v_str n-1 */
    char *inter_str = calloc(n_mult,1); 
    pointerEr(inter_str);
    
    i = 0;
    while (i < mult-1)
    {
      strcat(inter_str,v_str);
      strcat(inter_str,"->");
      i++;
    }
    strcat(inter_str,v_str);/* => itner_str = 2->2*/
    
    char *i_str = calloc(strlen(str)-strlen(b_str)+2,1);
    pointerEr(i_str);
    snprintf(i_str,strlen(str)-strlen(b_str)+1,"%s",str);
    
    char *new_str = 0;
    unsigned n_new_str = (unsigned)(strlen(i_str)+strlen(inter_str)+strlen(r_str))+1;
    new_str = calloc(n_new_str,1);
    pointerEr(new_str);
    
    strcat(new_str,i_str);
    strcat(new_str,inter_str);
    strcat(new_str,r_str);
    
    free(v_str);
    free(mult_str);
    free(m_str);
    free(inter_str);
    free(i_str);
    free(str);
    
    str = new_str;
    m_str = regex_find("\\(x[[:digit:]]+\\)",str);
  }
  ret_str = str;
  
  return ret_str;
}

/* having parameter name, it returns a pointer to 
// the corresponding parameter.
// ->return value: lhs=rhs of a parameter, 0 if not found
*/
Parameter_T *get_parameter(const char *const par_name)
{
  int i;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
      return parameters_global[i];
    i++;
  }
  
  return 0;
}

/* having the parameter name,
// it returns the double value of parameter which has been filled in double format.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: double value of parameter.
*/
double get_parameter_double_format(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  double v = DBL_MAX;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      if (!parameters_global[i]->double_flg)
        abortEr_s("Flag of parameter '%s' has not been set correctly.\n"
                  ,par_name);
        
      v = parameters_global[i]->rv_double;
      f = FOUND;
      break;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter '%s' couldn't be found.\n",par_name,file,line);
  }

  return v;
}

/* having the parameter name,
// it returns the array parameter which has been filled in array format.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: array value of parameter.
*/
double *get_parameter_array_format(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  double *v = 0;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = parameters_global[i]->rv_array;
      f = FOUND;
      break;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter '%s' couldn't be found.\n",par_name,file,line);
  }

  return v;
}


/* having the parameter name,
// it returns the INTEGER value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: integer value of parameter.
*/
int get_parameter_value_I(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  int v = INT_MAX;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = atoi(parameters_global[i]->rv);
      f = FOUND;
      break;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter '%s' couldn't be found.\n",par_name,file,line);
  }

  return v;
}

/* having the parameter name,
// it returns the DOUBLE value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: double value of parameter.
*/
double get_parameter_value_D(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  double v = DBL_MAX;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      if (parameters_global[i]->double_flg)
      {
        v = parameters_global[i]->rv_double;
      }
      else
      {
        v = strtod(parameters_global[i]->rv,0);
        parameters_global[i]->double_flg = 1;/* now we know this is double */
        parameters_global[i]->rv_double  = v;
      }
      f = FOUND;
      break;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter '%s' couldn't be found.\n",par_name,file,line);
  }

  return v;
}

/* having the parameter name,
// it returns the STRING value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: string value of parameter.
*/
const char *get_parameter_value_S(const char *const par_name,const char *const file, const int line,const Flag_T flg)
{
  char *v = 0;
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      v = parameters_global[i]->rv;
      f = FOUND;
      break;
    }
    i++;
  }
  
  if (flg == FATAL && f != FOUND)
  {
    abort_error_string("Parameter '%s' couldn't be found.\n",par_name,file,line);
  }
      
  return v;
}

/* reading the input file and make all of parameters 
// with their default value to the data base 
*/
int make_parameters(const char *const path)
{
  read_input_file(path);
  
  /* printing parameters */
  //if (test_print(PRINT_PARAMETERS))
    //pr_parameters();
  
  return EXIT_SUCCESS;
}

/* ->return value: count the total number of iterations */
unsigned total_iterations_ip(void)
{
  unsigned max = (unsigned)PgetiEZ("total_iterations_ip");
  char *subs = 0;
  unsigned l = 0,i;
  
  if (max != INT_MAX)/* if "total_iterations_ip" par is set */
    return max;
  
  i   = 0;
  max = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (parameters_global[i]->iterative)
    {
      l = 0;
      subs = strstr(parameters_global[i]->rv_ip,"->");
      l++;
      while (subs)
      {
        subs += 2;/* move forward */
        if (subs[0]=='\0')/* if after -> is empty */
          abortEr_s("No value is specified after '->' in parameter %s.\n",parameters_global[i]->lv);
        subs = strstr(subs,"->");
        if (subs)
          l++;
      }
      if (l > max)
        max = l;
    }
    i++;
  }
  
  Pseti("total_iterations_ip",(int)max+1);
  
  return max+1;/* +1 since the last value also is counted */
}

/* ->return value: count the total number of iterative parameters */
unsigned total_iterative_parameters_ip(void)
{
  unsigned i;
  unsigned count = 0;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (parameters_global[i]->iterative)
      count++;
    i++;
  }
  
  return count;
}

/* updating the iterative parameter for iteration number iter */
void update_iterative_parameter_ip(const unsigned iter)
{
  unsigned i;
  
  i = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (parameters_global[i]->iterative)
    {
      parameters_global[i]->rv = 
        get_n_value_str_ip(parameters_global[i],iter);
      if (parameters_global[i]->double_flg)
      {
        parameters_global[i]->rv_double  = 
          strtod(parameters_global[i]->rv,0);
      }
    }
    i++;
  }
    
}

/* getting the position number n-th of a value in the par and then
// return that value in string format, it frees the previous value and
// allocate memory for the new value.
// if the n is larger than the total number of values in the string
// it returns the last value.
// ->return value: the n-th value in iterative parameter in string format. */
char *get_n_value_str_ip(const Parameter_T *const par,const unsigned n)
{
  if (!par->iterative)
    abortEr_s("The parameter %s is not an iterative parameter",par->lv);
 
  if (!par->rv_ip)
    abortEr_s("The parameter %s doesn't have an iterative value",par->lv);
    
  char *ret = 0;
  char *subs,*subs2;
  const char *const rv_ip = par->rv_ip;
  unsigned j = 0;
  unsigned len = 0;
     
  if (par->rv)
    free(par->rv);
  
  j = 0;
  subs = strstr(rv_ip,"->");
  if (n == 0)
  {
    len = (unsigned)(subs-rv_ip+1);/* value-> */
    assert(len != UINT_MAX);
    ret = calloc(len,1);
    pointerEr(ret);
    strncpy(ret,rv_ip,len-1);
    ret[len-1] = '\0';
  }
  else
  {
    while (subs)
    {
      j++;
      subs += 2;
      subs2 = strstr(subs,"->");
      
      if (!subs2)/* if it is the last piece */
      {
        ret = dup_s(subs);
        break;
      }
      else if (j == n)
      {
        len = (unsigned)(subs2-subs+1);/* value1->value2 */
        assert(len != UINT_MAX);
        ret = calloc(len,1);
        pointerEr(ret);
        strncpy(ret,subs,len-1);
        ret[len-1] = '\0';
        break;
      }
      subs = subs2;
    }
  }
  
  return ret;   
}

/* ->return value: the name of n-th iterative variable, counting starts from 0 */
char *par_name_ip(const unsigned n)
{
  char *ret = 0;
  unsigned count;
  unsigned i;
  
  i = 0;
  count = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (parameters_global[i]->iterative)
    {
      if (count == n)
        return parameters_global[i]->lv;
      else
        count++;
    }
    i++;
  }
  
  if (count < n)
    abortEr("The total number of iterative parameters is fewer.\n");
  
  return ret;
}

/* ->return value: the value of n-th iterative variable in str format, counting starts from 0 */
char *par_value_str_ip(const unsigned n)
{
  char *ret = 0;
  unsigned count;
  unsigned i;
  
  i = 0;
  count = 0;
  while (parameters_global != 0 && parameters_global[i] != 0)
  {
    if (parameters_global[i]->iterative)
    {
      if (count == n)
        return parameters_global[i]->rv;
      else
        count++;
    }
    i++;
  }
  
  if (count < n)
    abortEr("The total number of iterative parameters is fewer.\n");
  
  return ret;
}

/* set the parameter unless it has already been set */
void set_default_parameter(const char *const lhs,const char *const rhs)
{
  const char *v;
  Parameter_T *par;
  
  par = get_parameter(lhs);
  if (par == 0)
    add_parameter(lhs,rhs);
  else
  {
    v = PgetsEZ(lhs);
    if (v == 0)
      par->rv = dup_s(rhs);
    else if (v[0] == '\0' || strcmp_i(v,"default"))
    {
      free(par->rv);
      par->rv = dup_s(rhs);
    }
  }
}

/* adding 2 block of memory for parameter data base 
// and putting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_parameter(Parameter_T ***const mem)
{
  unsigned i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  pointerEr((*mem));
  
  (*mem)[i] = calloc(1,sizeof(*(*mem)[i]));
  pointerEr((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* given the parameter name, free the parameter data base from it 
// and shrink the data base and put the last parameter in place of
// the deleted parameter. */
void free_parameter(const char *const par_name)
{
  Parameter_T *last_par = 0;
  unsigned np,i;
  
  /* count total number of parameters */
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  
  if (np == 0)
    return;
    
  for (i = 0; i < np; ++i)
  {
    if (strcmp_i(parameters_global[i]->lv,par_name))
    {
      last_par = parameters_global[np-1];
      
      _free(parameters_global[i]->lv);
      _free(parameters_global[i]->rv);
      _free(parameters_global[i]->rv_ip);
      _free(parameters_global[i]->rv_array);
      free(parameters_global[i]);
      
      parameters_global[i] = last_par;
      
      parameters_global = 
        realloc(parameters_global,np*sizeof(*parameters_global));
      pointerEr(parameters_global);
      
      parameters_global[np-1] = 0;
      break;
    }
  }
  
}

/* given the parameter, free its content and itself. 
// Note: it won't affect the parameter data base. */
void free_given_parameter(Parameter_T *par)
{
  if (!par)
    return;
    
  _free(par->lv);
  _free(par->rv);
  _free(par->rv_ip);
  _free(par->rv_array);
  free(par);
}

/* free the whole parameter data base */
void free_parameter_db(void)
{
  unsigned np;
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
  
    _free(parameters_global[np]->lv);
    _free(parameters_global[np]->rv);
    _free(parameters_global[np]->rv_ip);
    _free(parameters_global[np]->rv_array);
    free(parameters_global[np]);
    np++;
  }
  
  _free(parameters_global);
  
  parameters_global = 0;
}

