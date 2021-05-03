/*
// Alireza Rashti
// May 2018
*/

#include "parameters.h"

/* add or update a string value parameter */
void update_parameter_string(const char *const lv0, const char *const rv)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    Free(par->rv);
    par->rv = dup_s(rv);
    par->double_flg = 0;
  }
  else/* since it does not exist */
    add_parameter(lv,rv);
  
  Free(lv);
}

/* add or update an integer value parameter */
void update_parameter_integer(const char *const lv0, const int rv)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  char str_rv[100];
  sprintf(str_rv,"%d",rv);
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    Free(par->rv);
    par->rv = dup_s(str_rv);
    par->double_flg = 0;
  }
  else/* since it does not exist */
  {
    add_parameter(lv,str_rv);
  }
  
  Free(lv);
}

/* add or update a double value parameter, 
// if print_flg != 0 it prints the changes. */
void update_parameter_double_format(const char *const lv0, const double rv,const int print_flg)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  char str_rv[100] = {'\0'};
  
  par = get_parameter(lv);
  if (par)/* if this parameter exists update it */
  {
    if (print_flg)
    {
      char pr_msg[STR_SIZE1];
      double diff_a = 0.;/* absolute */
      double diff_r = 0.;/* relative */
      
      if (!EQL(par->rv_double,0.))
      {
        diff_a = (rv-par->rv_double);
        diff_r = diff_a/fabs(par->rv_double)*100.;
      }
        
      sprintf(pr_msg,PAR_FORMAT_PR,PAR_WIDTH_PR,lv,rv,diff_a,diff_r);
      printf(Pretty0"%s\n",pr_msg);
    }
    sprintf(str_rv,"%.20e",rv);
    Free(par->rv);
    
    /* NOTE:crucial to write in str format for checkpoint file purposes */
    par->rv         = dup_s(str_rv);
    par->rv_double  = rv;
    par->double_flg = 1;
  }
  else/* since it does not exist */
  {
    add_parameter_double(lv,rv,print_flg);
  }
  
  Free(lv);
}

/* adding left value and right value to parameter data base 
// double format. */
void add_parameter_double(const char *const lv0, const double rv,const int print_flg)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  char str_rv[100] = {'\0'};
  
  par = get_parameter(lv);
  if (par)
    Errors("This parameter \"%s\" has already been added!\n",lv);
    
  sprintf(str_rv,"%.20e",rv);
  par = alloc_parameter(&parameters_global);
  par->rv = dup_s(str_rv);
  par->lv = dup_s(lv);
  par->rv_double  = rv;
  par->double_flg = 1;
  
  if (print_flg)
  {
    char pr_msg[STR_SIZE1];
    const double diff = 0.;
    
    sprintf(pr_msg,PAR_FORMAT_PR,PAR_WIDTH_PR,lv,rv,diff,diff);
    printf(Pretty0"%s\n",pr_msg);
  }
  
  Free(lv);
}

/* adding left value and right value to parameter data base 
// array format. */
void add_parameter_array(const char *const lv0, const double *const rv,const Uint n)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  Uint i;
  
  par = get_parameter(lv);
  if (par)
    Errors("This parameter '%s' has already been added!\n",lv);
    
  par = alloc_parameter(&parameters_global);
  par->lv = dup_s(lv);
  
  par->rv_array = alloc_double(n);
  par->rv_n = n;
  for (i = 0; i < n; ++i)
    par->rv_array[i] = rv[i];
  
  Free(lv);
}

/* update parameter array format. */
void update_parameter_array(const char *const lv0, const double *const rv,const Uint n)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  Uint i;
  
  par = get_parameter(lv);
  
  if (par)
  {
    if (par->iterative)
      Errors("Wrong update: parameter '%s' is iterative.\n",lv);
    
    par->double_flg = 0;
      
    Free(par->rv_array);
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

  Free(lv);
}

/* adding left value and right value to parameter data base 
// string format */
void add_parameter_string(const char *const lv, const char *const rv)
{
  add_parameter(lv,rv);
}

/* adding left value and right value to parameter data base 
// string format (the most common one) */
void add_parameter(const char *const lv0, const char *const rv)
{
  Parameter_T *par;
  char *lv = PAR_NAME_RULE(lv0);
  
  par = get_parameter(lv);
  if (par)
    Errors("This parameter \"%s\" has already been added!\n",lv);
    
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
      if (par->rv[0] == 0)
        Errors("No value for parameter '%s' .\n",par->lv);

    }
    else
    {
      par->iterative = 0;
      par->rv = dup_s(rv);
    }
  }
  
  Free(lv);
}

/* some iterative parameter may have multiplicity with (x?), e.g:
// par = 1->2(x2)->5(x3)->8
// this function parses it and returns it in regular format, i.e:
// par = 1->2->2->5->5->5->8
// ->return value: an iterative parameter with regular format */
static char *parse_multiplicity_of_iterative_parameter(const char *const rv)
{
  if (regex_search("->$",rv))
    Errors("Wrong syntax for '%s'; '->' at the end of the line.",rv);
    
  char *ret_str = 0;
  char *str = dup_s(rv);
  
  char *m_str = regex_find("\\(x[[:digit:]]+\\)",rv);
  while (m_str)
  {
    /* reference: 1->2(x2)->5(x3)->8 */
    Uint len = (Uint)strlen(m_str);/* (x2) => 4 */
    const char *l_str = strstr(str,m_str);/*    => (x2)->5(x3)->8 */
    const char *r_str = strlen(l_str) > len ? l_str+len: "\0";/* => ->5(x3)->8 */
    
    const char *b_str = l_str;
    while (b_str != str && *b_str != '>')
      b_str--;
    /* b_str => >2 */  
    if (*b_str == '>')
      b_str++;/* b_str => 2 */
    
    Uint n = (Uint)(l_str-b_str+1);  
    char *v_str = calloc(n,1);
    IsNull(v_str);
    
    Uint i = 0;
    while (i < n-1)
    {
      v_str[i] = b_str[i];
      i++;
    }/* =? v_str = 2 */
    v_str[n-1] = '\0';
    
    char *mult_str  = regex_find("[[:digit:]]+",m_str);
    Uint mult   = (Uint)atoi(mult_str);
    
    if (mult == UINT_MAX)
      Errors("Wrong syntax for '%s'; negative multiplicity.\n",rv);
    if (mult == 0)
      Errors("Wrong syntax for '%s'; zero multiplicity.\n",rv);
    
    Uint n_mult = (mult-1)*2+mult*(n-1)+1;/* for each '->' 2 Byte,
                                    // for each v_str n-1 */
    char *inter_str = calloc(n_mult,1); 
    IsNull(inter_str);
    
    i = 0;
    while (i < mult-1)
    {
      strcat(inter_str,v_str);
      strcat(inter_str,"->");
      i++;
    }
    strcat(inter_str,v_str);/* => itner_str = 2->2*/
    
    char *i_str = calloc(strlen(str)-strlen(b_str)+2,1);
    IsNull(i_str);
    snprintf(i_str,strlen(str)-strlen(b_str)+1,"%s",str);
    
    char *new_str = 0;
    Uint n_new_str = (Uint)(strlen(i_str)+strlen(inter_str)+strlen(r_str))+1;
    new_str = calloc(n_new_str,1);
    IsNull(new_str);
    
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
Parameter_T *get_parameter(const char *const par_name0)
{
  if (!parameters_global)
    return 0;
    
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
    {
      Free(par_name);
      return parameters_global[i];
    }
    i++;
  }
  
  Free(par_name);
  
  return 0;
}

/* having the parameter name,
// it returns the double value of parameter which has been filled in double format.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: double value of parameter.
*/
double get_parameter_double_format(const char *const par_name0,const char *const file, const int line,const Flag_T flg)
{
  double v = DBL_MAX;
  
  if (!parameters_global)
    return v;
  
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
    {
      if (!parameters_global[i]->double_flg)
        Errors("Flag of parameter '%s' has not been set correctly.\n"
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

  Free(par_name);
  
  return v;
}

/* having the parameter name,
// it returns the array parameter which has been filled in array format.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: array value of parameter.
*/
double *get_parameter_array_format(const char *const par_name0,const char *const file, const int line,const Flag_T flg)
{
  double *v = 0;
  
  if (!parameters_global)
    return v;
  
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
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

  Free(par_name);
  
  return v;
}


/* having the parameter name,
// it returns the INTEGER value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: integer value of parameter.
*/
int get_parameter_value_I(const char *const par_name0,const char *const file, const int line,const Flag_T flg)
{
  int v = INT_MAX;
  
  if (!parameters_global)
    return v;
  
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
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

  Free(par_name);
  return v;
}

/* having the parameter name,
// it returns the DOUBLE value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: double value of parameter.
*/
double get_parameter_value_D(const char *const par_name0,const char *const file, const int line,const Flag_T flg)
{
  double v = DBL_MAX;
  
  if (!parameters_global)
    return v;
  
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
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

  Free(par_name);
  return v;
}

/* having the parameter name,
// it returns the STRING value of parameter.
// if flag == FATAL and couldn't find the par_name, gives error.
// ->return value: string value of parameter.
*/
const char *get_parameter_value_S(const char *const par_name0,const char *const file, const int line,const Flag_T flg)
{
  char *v = 0;
  
  if (!parameters_global)
    return v;
  
  char *par_name = PAR_NAME_RULE(par_name0);
  int i;
  Flag_T f = NONE;
  
  i = 0;
  while (parameters_global[i])
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
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

  Free(par_name);
  
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
Uint total_iterations_ip(void)
{
  Uint max = (Uint)PgetiEZ("total_iterations_ip");
  char *subs = 0;
  Uint l = 0,i;
  
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
          Errors("No value is specified after '->' in parameter %s.\n",parameters_global[i]->lv);
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
Uint total_iterative_parameters_ip(void)
{
  Uint i;
  Uint count = 0;
  
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
void update_iterative_parameter_ip(const Uint iter)
{
  Uint i;
  
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
char *get_n_value_str_ip(const Parameter_T *const par,const Uint n)
{
  if (!par->iterative)
    Errors("The parameter %s is not an iterative parameter",par->lv);
 
  if (!par->rv_ip)
    Errors("The parameter %s doesn't have an iterative value",par->lv);
    
  char *ret = 0;
  char *subs,*subs2;
  const char *const rv_ip = par->rv_ip;
  Uint j = 0;
  Uint len = 0;
     
  if (par->rv)
    free(par->rv);
  
  j = 0;
  subs = strstr(rv_ip,"->");
  if (n == 0)
  {
    len = (Uint)(subs-rv_ip+1);/* value-> */
    assert(len != UINT_MAX);
    ret = calloc(len,1);
    IsNull(ret);
    strncpy(ret,rv_ip,len-1);
    ret[len-1] = '\0';
  }
  else
  {
    if (subs)
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
        len = (Uint)(subs2-subs+1);/* value1->value2 */
        assert(len != UINT_MAX);
        ret = calloc(len,1);
        IsNull(ret);
        strncpy(ret,subs,len-1);
        ret[len-1] = '\0';
        break;
      }
      subs = subs2;
    }
    else/* it might happen that user written (x1) thus no '->' exists */
    {
      ret = dup_s(rv_ip);
    }
  }
  
  return ret;   
}

/* ->return value: the name of n-th iterative variable, counting starts from 0 */
char *par_name_ip(const Uint n)
{
  char *ret = 0;
  Uint count;
  Uint i;
  
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
    Error0("The total number of iterative parameters is fewer.\n");
  
  return ret;
}

/* ->return value: the value of n-th iterative variable in str format, counting starts from 0 */
char *par_value_str_ip(const Uint n)
{
  char *ret = 0;
  Uint count;
  Uint i;
  
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
    Error0("The total number of iterative parameters is fewer.\n");
  
  return ret;
}

/* set the parameter unless it has already been set */
void set_default_parameter(const char *const lhs0,const char *const rhs)
{
  const char *v;
  Parameter_T *par;
  char *lhs = PAR_NAME_RULE(lhs0);

  par = get_parameter(lhs);
  if (par == 0)
    add_parameter(lhs,rhs);
  else
  {
    v = PgetsEZ(lhs);
    if (v == 0)
      par->rv = dup_s(rhs);
    else if (v[0] == '\0')
    {
      free(par->rv);
      par->rv = dup_s(rhs);
    }
  }
  
  Free(lhs);
}

/* adding 2 block of memory for parameter data base 
// and putting the last block to null and 
// returning pointer to one before the last block
*/
void *alloc_parameter(Parameter_T ***const mem)
{
  Uint i;
  
  for (i = 0; (*mem) != 0 && (*mem)[i] != 0 ; i++);
  
  (*mem) = realloc((*mem),(i+2)*sizeof(*(*mem)));
  IsNull((*mem));
  
  (*mem)[i] = calloc(1,sizeof(*(*mem)[i]));
  IsNull((*mem)[i]);
  
  (*mem)[i+1] = 0;
  
  return (*mem)[i];
}

/* given the parameter name, free the parameter data base from it 
// and shrink the data base and put the last parameter in place of
// the deleted parameter. */
void free_parameter(const char *const par_name0)
{
  Parameter_T *last_par = 0;
  char *par_name = PAR_NAME_RULE(par_name0);
  Uint np,i;
  
  /* count total number of parameters */
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
    np++;
  
  if (np == 0)
  {
    Free(par_name);
    return;
  }
    
  for (i = 0; i < np; ++i)
  {
    if (!strcmp(parameters_global[i]->lv,par_name))
    {
      last_par = parameters_global[np-1];
      
      Free(parameters_global[i]->lv);
      Free(parameters_global[i]->rv);
      Free(parameters_global[i]->rv_ip);
      Free(parameters_global[i]->rv_array);
      free(parameters_global[i]);
      
      parameters_global[i] = last_par;
      
      parameters_global = 
        realloc(parameters_global,np*sizeof(*parameters_global));
      IsNull(parameters_global);
      
      parameters_global[np-1] = 0;
      break;
    }
  }
  
  Free(par_name);
}

/* given the parameter, free its content and itself. 
// Note: it won't affect the parameter data base. */
void free_given_parameter(Parameter_T *par)
{
  if (!par)
    return;
    
  Free(par->lv);
  Free(par->rv);
  Free(par->rv_ip);
  Free(par->rv_array);
  free(par);
}

/* free the whole parameter data base */
void free_parameter_db(void)
{
  Uint np;
  
  np = 0;
  while (parameters_global != 0 && parameters_global[np] != 0)
  {
  
    Free(parameters_global[np]->lv);
    Free(parameters_global[np]->rv);
    Free(parameters_global[np]->rv_ip);
    Free(parameters_global[np]->rv_array);
    free(parameters_global[np]);
    np++;
  }
  
  Free(parameters_global);
  
  parameters_global = 0;
}

/* ->: if exceeds max iteration returns 1, otherwise 0.
// updating iterative parameters, STOP parameter and output directories.
// new output directory is made based on changing of resolution. */
int 
update_iteration_params
  (
   const Uint main_loop_iter,
   const char *const prefix/* parameter prefix */,
   const char *const dir_name_format/* eg: "BHNS_%s_%ux%ux%u" */
  )
{
  assert(dir_name_format);
  char prepar[PAR_LEN] = {'\0'};
  const Uint total_iters = total_iterations_ip();
  const Uint total_ipars = total_iterative_parameters_ip();
  const Uint iter_n     = (Uint)Pgeti(PrefixIt(prefix,"iteration"));
  const Uint res_iter_n = (Uint)Pgeti(PrefixIt(prefix,"resolution_iteration"));
  Uint iter;/* number of iterations have been performed for the simulation */
  Uint n[3];/* number of points */
  char folder_name_next[STR_SIZE2] = {'\0'},
       folder_name_prev[STR_SIZE2] = {'\0'};
  char *folder_path,*folder_path2;
  const char *const parfile_name = Pgets("parameter_file_name");/* no prefix */
  const char *const parfile_stem = Pgets("parameter_file_name_stem");/* no prefix */
  char cp_cmd[STR_SIZE2];
  char *str;
  Uint i;
  
  /* if the start is from checkpoint_file do nothing */
  //if (Pcmps(PrefixIt(prefix,"start_off"),"checkpoint_file"))
  //{
    //Pset_default("top_directory","NOT_SPECIFIED_YET");
    //return 0;
  //}
  
  /* if top directory is not set */
  if (!PgetsEZ("top_directory"))
  {
    folder_path = make_directory(Pgets("relative_root_path"),parfile_stem);
    Pset_default("top_directory",folder_path);
    Free(folder_path);
  }
  
  /* when starting from checkpoint, iter_n > main_loop_iter 
  // so to avoid redo the simulations we set iter to the largest. */
  if (iter_n > main_loop_iter)
  {
    iter = iter_n;
    /* updating iterative parameters for the new round of iteration */
    Pseti(PrefixIt(prefix,"iteration"),(int)iter+1);/* +1 is crucial */
  }
  else
  {
    iter = main_loop_iter;
    /* updating iterative parameters for the new round of iteration */
    Pseti(PrefixIt(prefix,"iteration"),(int)iter);
  }
  
  /* if exceeds total iteration => stop */
  if (iter >= total_iters)
  {
    Pseti(PrefixIt(prefix,"STOP"),1);
    Pseti(PrefixIt(prefix,"iteration"),(int)iter);
    return 1;
  } 
  
  /* find the previous folder name */
  n[0] = (Uint)PgetiEZ("n_a");
  n[1] = (Uint)PgetiEZ("n_b");
  n[2] = (Uint)PgetiEZ("n_c");
  sprintf(folder_name_prev,dir_name_format,parfile_stem,n[0],n[1],n[2]);  
  
  /* update the parameter accoding to the iteration number */
  update_iterative_parameter_ip(iter);
  
  /* print the iterative parameters */
  if (total_ipars && main_loop_iter)
  {
    printf("Next iterative parameter(s) are:\n");
    for (i = 0; i < total_ipars; ++i)
      printf(Pretty0"%-30s = %-15s\n",par_name_ip(i),par_value_str_ip(i));
  }
  
  /* find the name of next folder */
  n[0] = (Uint)PgetiEZ("n_a");
  n[1] = (Uint)PgetiEZ("n_b");
  n[2] = (Uint)PgetiEZ("n_c");
  
  /* this parameter helps to use some of the previous grid data */
  Pseti(PrefixIt(prefix,"did_resolution_change?"),0);
  
  /* update resolution iter. */
  Pseti(PrefixIt(prefix,"resolution_iteration"),(int)(res_iter_n+1));
  
  sprintf(folder_name_next,dir_name_format,parfile_stem,n[0],n[1],n[2]);
  /* if the resolution isn't the same or it is the first iteration */
  if (strcmp(folder_name_next,folder_name_prev) || iter == 0)/* if n is updated */
  {
    /* iteration number used in solving, reset this for each resolution */
    Pseti(PrefixIt(prefix,"resolution_iteration"),0);
    sprintf(folder_name_next,dir_name_format,parfile_stem,n[0],n[1],n[2]);
    folder_path = make_directory(Pgets("top_directory"),folder_name_next);
    Psets(PrefixIt(prefix,"my_directory"),folder_path);
    folder_path2 = make_directory(folder_path,"Diagnostics");
    Psets(PrefixIt(prefix,"Diagnostics"),folder_path2);
    
    /* copy the parameter file into the new directory */
    str = strrchr(folder_path,'/');
    assert(str);
    str++;
    sprintf(cp_cmd,"cp %s/%s %s/%s.par",
    Pgets("relative_root_path")/* no prefix */,
    parfile_name,folder_path,str);
    shell_command(cp_cmd);
    
    free(folder_path);
    free(folder_path2);
    
    /* => resolution changed */
    Pseti(PrefixIt(prefix,"did_resolution_change?"),1);
  }
  
  return 0;
}

/* ->: alloc memory and return ruled1 s, if s==0 returns 0.
// the rule(s) for parameter name regarding lower case or upper case.
// thus, there would be no convert to lower case to get parameter.
//
// rule1:
// ======
// if there is '_' in string s, then the first piece of str 
// from begining to the first '_' become upper case 
// and the rest lower case.
// if there is no '_' in s, make it all lower case. */
static char *par_rule1_uppercase_lowercase(const char *const s)
{
  assert(s);
  
  char *rule_s = dup_s(s);/* note: the last char is '\0' */
  char *aux = 0;
  Uint i    = 0;
 
  aux = strchr(rule_s,'_');
  if (aux)
  {
    /* make it upper case till '_' */
    i = 0;
    while(&rule_s[i] != aux && 
           rule_s[i] != '\0')
    {
      rule_s[i] = (char)toupper(rule_s[i]);
      i++;
    }
    /* the rest is lower case */
    while(rule_s[i] != '\0')
    {
      rule_s[i] = (char)tolower(rule_s[i]);
      i++;
    }
  }
  else
  {
    i = 0;
    while(rule_s[i] != '\0')
    {
      rule_s[i] = (char)tolower(rule_s[i]);
      i++;
    }
  }
  
  return rule_s;
}

/* ->: alloc memory and return ruled2 string s, if s == 0 returns 0.
// the rule(s) for parameter name regarding lower case or upper case.
// thus, there would be no convert to lower case to get parameter.
//
// rule2:
// ======
// convert the whole string s to lowercase. */
static char *par_rule2_lowercase(const char *const s)
{
  assert(s);
  
  char *const rule_s = dup_s(s);/* note: the last char is '\0'. */
  Uint i = 0;
  
  i = 0;
  while(rule_s[i] != '\0')
  {
    rule_s[i] = (char)tolower(rule_s[i]);
    i++;
  }
  
  return rule_s;
}

/* ->: casted parameter's name according to the rule */
char *par_name_rule(const char *const s)
{
  return PAR_NAME_RULE(s);
}

