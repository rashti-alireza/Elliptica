/*
// Alireza Rashti
// July 2018
*/

/* Print Fields:
// note: in Silo language, I'm using curvilinear format for both mesh and data (fields)
// regardless of the patch is Cartesian or not.
// note: data and mesh in Silo must be written in column major order,
//       otherwise for inhomogeneous resolutions you get a scrambled mesh.
//
// usage examples:
// ===============
// # parameter that is determined in input file is like:
// silo_output_3d = (V_U0,V_U1,V_U2),psi,eta,(a_U0,a_U1,a_U2)
// # as one can see the vector quantities determined by parenthesis 
//
// Pr_Field_T *pr  = init_PrField(grid);
// pr->folder      = "folder_path";
// pr->cycle       = iteration_number;// if you wanna plot data at each iteration
//
// # the following options and flag are not necessary, their default value is 0.
// pr->multimesh_f = 1; # if you wanna make a master file for all patches as a whole grid
// pr->multivar_f  = 1; # if you wanna make a master file for all fields. This, option seems not working with VisIt.
// pr->abc_f       = 1; # if you wanna have the patches and fields in (X,Y,Z) corods or (a,b,c) coords.
//
// # print all patchs and fields
// pr_fields(pr);
//
// # free
// free_PrField(pr);
*/

#include "pr_for_fields.h"

/* DeLimits */
#define DL_OP '('
#define DL_CP ')'
#define DL_C ','

/* given print parameter related to fields, the folder, and time = cycle,
// it reads the parameter and the fields indicated there 
// and prints the whole grid and fields in the specified folder. */
void pr_fields(Pr_Field_T *const pr)
{
  if (!pr)
    return;
  
  /* using silo */  
  if (PgetsEZ("silo_output_3d"))
  {
    parse_parameter_3d(Pgets("silo_output_3d"),pr);
    pr_hdf5_silo(pr);
  }
  
}
/* initiating a Pr_Field_T for printing with given grid.
// ->return value: Pr_Field_T */
Pr_Field_T *init_PrField(const Grid_T *const grid)
{
  Pr_Field_T *pr = calloc(1,sizeof(*pr));
  IsNull(pr);
  pr->grid = grid;
  
  return pr;
}

/* freeing Pr_Field_T */
void free_PrField(Pr_Field_T *pr)
{
  struct Info_S *info = pr->group;
  Uint i;
  
  for (i = 0; i < pr->ng; ++i)
  {
    if (info[i].vec_flg)
    {
      free(info[i].comp[0]);
      free(info[i].comp[1]);
      free(info[i].comp[2]);
    }
    else
      free(info[i].field);
    
  }
  Free(info);
  Free(pr);
}

/* it reads the parameter to find out what and how it's going to be print.
// given parameter it finds the name of all of the fields to be printed 
// and the coordinates which these field evaluated on.
// coords can be Cartesian x,y,z or Curvilinear a,b,c.
// NOTE: fields are supposed to be written like {(field1)vs(a,b,c)|coord}. 
// NOTE: it supports regular expression for scalar, ex: ^dpsi_D.+ */
static void parse_parameter_3d(const char *const par,Pr_Field_T *const pr)
{
  struct Info_S *info_s = 0,*Pinfo;
  char *tok = dup_s(par);
  char *save = 0,*sub_tok = 0;
  Uint Ninfo = 0;
  
  save    = tok;
  sub_tok = tok_s(sub_tok,DL_C,&save);
  /* either => sub_tok = f1 and save = f2,f3,... 
  // or     => sub_tok = (fx,fy,fz) and save = ... */
  if (sub_tok == 0)
    fprintf(stderr,"No field in parameter file is specified for 3d printing.\n");
    
  while (sub_tok)
  {
    char *savess = 0,*ss = 0,*dump;
    
    info_s = realloc(info_s,(Ninfo+1)*sizeof(*info_s)); 
    IsNull(info_s);
    Pinfo = &info_s[Ninfo];
    Pinfo->field   = 0;
    Pinfo->comp[0] = 0;
    Pinfo->comp[1] = 0;
    Pinfo->comp[2] = 0;
    
    /* if sub_tok = (fx  save = fy,fz)... */
    if (strchr(sub_tok,DL_OP))
    {
      savess = save;
      save = strchr(save,DL_CP);
      save++;
      if (save[0] == DL_C)
        save++;
        
      ss = tok_s(sub_tok,DL_OP,&dump);/* => ss = fx,dump = '\0' */
      Pinfo->comp[0] = dup_s(ss);
      ss = tok_s(0,DL_C,&savess); /* => ss = fy, savess = fz) */
      Pinfo->comp[1] = dup_s(ss);
      ss = tok_s(0,DL_CP,&savess); /* => ss = fz */
      Pinfo->comp[2] = dup_s(ss);
      Pinfo->vec_flg = 1;
    }
    /* sub_tok = f1 */
    else
    {
      Pinfo->field = dup_s(sub_tok);
      Pinfo->vec_flg = 0;
    }
    sub_tok = tok_s(0,DL_C,&save);
    Ninfo++;
  }
  
  pr->group = info_s;
  pr->ng = Ninfo;
  
  free(tok);
}

/* a quick print of given grid and all of the fields determined
// in the parameter file.*/
int print_fields_3D(const Grid_T *const grid,const int iteration,
                    const char *const dir)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  Pr_Field_T *pr  = init_PrField(grid);
  pr->folder      = dir;
  pr->cycle       = iteration;
  pr_fields(pr);
  free_PrField(pr);

  FUNC_TOC
  return EXIT_SUCCESS;
}

/* ->: largest L2 norm.
// print norms L2, L1 and L_inf of the specified fields 
// in parameter "txt_output_0d" (supporting regular expression) */
double print_fields_0D(const Grid_T *const grid,const int iteration, 
                      const char *const folder)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  /* list of the fields to be printed out */
  char **f = 
     read_separated_items_in_string(PgetsEZ("txt_output_0d"),',');
  double largest_L2_error = 0; 
  Uint i,p;
  
  if(f)         
  for (i = 0; f[i]; ++i)
  {
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      Field_T *field = 0;
      Uint nn        = patch->nn;
      int field_ind  = _Ind(f[i]);
      Uint Nf        = 0;
      Uint *fInd     = 0;
      Uint Ufield_ind= 0;
      Uint fi;
      double Linf,L2,L1;
      FILE *file_Linf,*file_L1,*file_L2;
      char file_name_Linf[STR_LEN];
      char file_name_L1[STR_LEN];
      char file_name_L2[STR_LEN];
      char *stem = strstr(patch->name,"_");
      stem++;
      
      /* if couldn't find, maybe its a regex. */
      if (field_ind < 0)
      {
        fInd = find_field_index_regex(patch,f[i],&Nf);
        
        /* if nothing found. */
        if (!fInd)
          continue;
      }
      else
      {
        Ufield_ind = (Uint)field_ind;
        fInd = &Ufield_ind;
        Nf   = 1;
      }
      
      for (fi = 0; fi < Nf; ++fi)
      {
        field = patch->fields[fInd[fi]];
        if (!field->v)/* if field empty */
          continue;
          
        sprintf(file_name_Linf,"%s/%s_Linf_%s.txt",folder,field->name,stem);
        sprintf(file_name_L1,  "%s/%s_L1_%s.txt" ,folder,field->name,stem);
        sprintf(file_name_L2,  "%s/%s_L2_%s.txt" ,folder,field->name,stem);
        
        if (access(file_name_Linf,F_OK) != -1)/* if file exists */
        {
          file_Linf = Fopen(file_name_Linf,"a");
        }
        else
        {
          file_Linf = Fopen(file_name_Linf,"w");
          fprintf(file_Linf,"#iteration  %s\n",field->name);
        }
        
        if (access(file_name_L1,F_OK) != -1)
        {
          file_L1 = Fopen(file_name_L1,"a");
        }
        else
        {
          file_L1 = Fopen(file_name_L1,"w");
          fprintf(file_L1,"#iteration  %s\n",field->name);
        }
        
        if (access(file_name_L2,F_OK) != -1)
        {
          file_L2 = Fopen(file_name_L2,"a");
        }
        else
        {
          file_L2 = Fopen(file_name_L2,"w");
          fprintf(file_L2,"#iteration  %s\n",field->name);
        }
          
        Linf  = L_inf(nn,field->v);
        L2    = L2_norm(nn,field->v,0);
        L1    = L1_norm(nn,field->v,0);
        
        fprintf(file_Linf,"%-11u %0.15f\n",iteration,Linf);
        fprintf(file_L1,  "%-11u %0.15f\n",iteration,L1);
        fprintf(file_L2,  "%-11u %0.15f\n",iteration,L2);
        
        largest_L2_error = L2 > largest_L2_error ? L2 : largest_L2_error;
        
        Fclose(file_Linf);
        Fclose(file_L1);
        Fclose(file_L2);
      }
      /* only if it is a regex found free it. */
      if (field_ind < 0)
        Free(fInd);
    }
  }
  free_2d(f);
  
  FUNC_TOC
  return largest_L2_error;
}

/* print the value of the given fields in the params along a 
// specified line.
// note: parameter "txt_output_1d" supports regular expression. */
void print_fields_1D(const Grid_T *const grid,const int iteration, 
                      const char *const folder)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  /* parsed directions. */
  struct parsed_s
  {
    Uint Xline:1;/* if X-line is asked 1, otherwise 0. */
    Uint Yline:1;/* if Y-line is asked 1, otherwise 0. */
    Uint Zline:1;/* if Z-line is asked 1, otherwise 0. */
    double X;/* X coords of the print line */
    double Y;/* Y coords of the print line */
    double Z;/* Z coords of the print line */
    char strv[STR_LEN];/* save a version of param values. */
    char suffix[STR_LEN/2];/* suffix for the file. */
  };
  const char *const line_par_name = "txt_output_1d_line";
  /* list of the fields(flds) to be printed out */
  char **flds = 
     read_separated_items_in_string(PgetsEZ("txt_output_1d"),',');
  /* list of requested directions(d) */
  char **d = 
     read_separated_items_in_string(PgetsEZ(line_par_name),')');
  
  Uint Nparsed;/* number of parsed struct */
  Uint i,j,k,f,l,p,counter;
  
  if (!d)
    printf(Pretty0"No line was specified.\n");
  
  /* count the number of meaning full values */
  counter = 0;
  for (i = 0; d && d[i]; ++i)
  {
    if (strlen(d[i]) > 5)/* coz we need, for example, "(I,1,2". 
                         // note ')' is eliminated. */
      counter++;
  }
  Nparsed = counter;
  struct parsed_s parsed[counter];
  j = 0;
  for (i = 0; d && d[i]; ++i)
  {
     counter = 0;
     if (strlen(d[i]) <= 5)
       continue;
     /* substitue ,e.g., ",(0,J,2" by "0,J,2". */
     regex_replace(d[i],"^.?\\(","",d[i]);
     char **subs = read_separated_items_in_string(d[i],',');
     
     /* sanity checks */
     if (!subs) Errors("Wrong format for %s.\n",line_par_name);
     for (counter = 0; subs[counter]; counter++);
     if (counter != 3) Errors("Wrong format for %s.\n",line_par_name);
     
     if (
         /* 'I/i' should only take place on 0th position */
         *subs[1] == 'X' || 
         *subs[1] == 'x' || 
         *subs[2] == 'X' || 
         *subs[2] == 'x' ||
         /* 'J/j' should only take place on 1st position */
         *subs[0] == 'Y' || 
         *subs[0] == 'y' || 
         *subs[2] == 'Y' || 
         *subs[2] == 'y' || 
         /* 'K/k' should only take place on 2nd position */
         *subs[0] == 'Z' || 
         *subs[0] == 'z' || 
         *subs[1] == 'Z' || 
         *subs[1] == 'z'
         )
     {
       Errors("Wrong order for %s.\n",line_par_name); 
     }
     
     /* realize */
     counter = 0;
     if (*subs[0] == 'X' || *subs[0] == 'x')
     {
       parsed[j].Xline = 1;
       parsed[j].Yline = 0;
       parsed[j].Zline = 0;
       
       parsed[j].X = DBL_MAX;
       parsed[j].Y = atof(subs[1]);
       parsed[j].Z = atof(subs[2]);
       sprintf(parsed[j].suffix,"X_%s_%s",subs[1],subs[2]);
       counter++;
     }
     if (*subs[1] == 'Y' || *subs[1] == 'y')
     {
       parsed[j].Xline = 0;
       parsed[j].Yline = 1;
       parsed[j].Zline = 0;
       
       parsed[j].X = atof(subs[0]);
       parsed[j].Y = DBL_MAX;
       parsed[j].Z = atof(subs[2]);
       sprintf(parsed[j].suffix,"%s_Y_%s",subs[0],subs[2]);
       counter++;
     }
     if (*subs[2] == 'Z' || *subs[2] == 'z')
     {
       parsed[j].Xline = 0;
       parsed[j].Yline = 0;
       parsed[j].Zline = 1;
       
       parsed[j].X = atof(subs[0]);
       parsed[j].Y = atof(subs[1]);
       parsed[j].Z = DBL_MAX;
       sprintf(parsed[j].suffix,"%s_%s_Z",subs[0],subs[1]);
       counter++;
     }
     /* if we have more than two letters, e.g., "(I,J,7)". */
     if (counter != 1)
       Errors("Wrong order for %s.\n",line_par_name); 
       
     sprintf(parsed[j].strv,"(%s)",d[i]);
     j++;
     free_2d(subs);
  }
  
  /* test */
  if (1)
  for (i = 0; i < Nparsed; ++i)
  {
    printf("param[%u] = %s => (%u,%u,%u) & (%g,%g,%g)\n",
            i, parsed[i].strv, 
            parsed[i].Xline,parsed[i].Yline,parsed[i].Zline,
            parsed[i].X,parsed[i].Y,parsed[i].Z);
  }
  
  if(flds && d)
  for (l = 0; l < Nparsed; ++l)
  {
    struct parsed_s *line = &parsed[l];
    Uint ijk;/* (X,Y,Z) indices */
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      Field_T **fields = 0;/* save the found fields */
      Uint Nfld = 0;/* total num of the found fields. */
      Field_T *field = 0;
      //Uint nn        = patch->nn;
      int field_ind  = 0;
      Uint Nf        = 0;
      Uint *fInd     = 0;
      Uint Ufield_ind= 0;
      Uint fi;
      
      const Uint *const n = patch->n;
      Uint I,J,K;
      double min_d;
      FILE *file;
      char file_name[STR_LEN];
      char *stem = strstr(patch->name,"_");
      stem++;
      
      /* finding the index of the closest point to 
      // the specified one in the param. */
      I = J = K = UINT_MAX;
      if (line->Xline)
      {
        min_d = DBL_MAX;
        for (j = 0; j < n[1]; ++j)
        {
          ijk = i_j_k_to_ijk(n,0,j,0);
          double dist = fabs(line->Y-patch->node[ijk]->X[1]);
          if (dist < min_d)
          {
            J = j;
            min_d = dist;
          }
        }
        
        min_d = DBL_MAX;
        for (k = 0; k < n[2]; ++k)
        {
          ijk = i_j_k_to_ijk(n,0,0,k);
          double dist = fabs(line->Z-patch->node[ijk]->X[2]);
          if (dist < min_d)
          {
            K = k;
            min_d = dist;
          }
        }
      }
      
      else if (line->Yline)
      {
        min_d = DBL_MAX;
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,0,0);
          double dist = fabs(line->X-patch->node[ijk]->X[0]);
          if (dist < min_d)
          {
            I = i;
            min_d = dist;
          }
        }
        
        min_d = DBL_MAX;
        for (k = 0; k < n[2]; ++k)
        {
          ijk = i_j_k_to_ijk(n,0,0,k);
          double dist = fabs(line->Z-patch->node[ijk]->X[2]);
          if (dist < min_d)
          {
            K = k;
            min_d = dist;
          }
        }
      }
      
      else if (line->Zline)
      {
        min_d = DBL_MAX;
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,0,0);
          double dist = fabs(line->X-patch->node[ijk]->X[0]);
          if (dist < min_d)
          {
            I = i;
            min_d = dist;
          }
        }
        
        min_d = DBL_MAX;
        for (j = 0; j < n[1]; ++j)
        {
          ijk = i_j_k_to_ijk(n,0,j,0);
          double dist = fabs(line->Y-patch->node[ijk]->X[1]);
          if (dist < min_d)
          {
            J = j;
            min_d = dist;
          }
        }
      }
      else
      {
        Error0(NO_OPTION);
      }
      
      /* find all fields for this patch */
      for (f = 0; flds[f]; ++f)
      {
        Nf   = 0;
        fInd = 0;
        Ufield_ind = 0;
        field_ind  = _Ind(flds[f]); 
        /* if couldn't find, maybe its a regex. */
        if (field_ind < 0)
        {
          fInd = find_field_index_regex(patch,flds[f],&Nf);
            
          /* if nothing found. */
          if (!fInd)
            continue;
        }
        else
        {
          Ufield_ind = (Uint)field_ind;
          fInd = &Ufield_ind;
          Nf   = 1;
        }
        for (fi = 0; fi < Nf; ++fi)
        {
          field = patch->fields[fInd[fi]];
          if (!field->v)/* if field empty */
            continue;
          
          fields = realloc(fields,(Nfld+1)*sizeof(*fields));
          IsNull(fields);
          fields[Nfld] = field;
          Nfld++;
        }
        /* only if it is a regex found free it. */
        if (field_ind < 0)
        Free(fInd);
      }/* end of for (f = 0; flds[f]; ++f) */
      
      /* create the file if not exist */
      sprintf(file_name,"%s/%s_%s_1d.txt",folder,stem,line->suffix);
      if (access(file_name,F_OK) != -1)/* if file exists */
      {
        file = Fopen(file_name,"a");
      }
      /* open a new file with an appropriate header */
      else
      {
        file = Fopen(file_name,"w");
        
        /* header:  iteration field1 field2 ... */
        fprintf(file,"# iteration");
        for (f = 0; f < Nfld; ++f)
          fprintf(file," %s", field->name);
      }
      fprintf(file,"# \"time = %d\" %s\n",iteration,line->strv);
      
      /* print coords and fields in each column. */
      if (line->Xline)
      {
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,J,K);
          fprintf(file," %0.15f ",patch->node[ijk]->X[0]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      else if (line->Yline)
      {
        for (j = 0; j < n[1]; ++j)
        {
          ijk = i_j_k_to_ijk(n,I,j,K);
          fprintf(file," %0.15f ",patch->node[ijk]->X[1]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      else if (line->Zline)
      {
        for (k = 0; k < n[2]; ++k)
        {
          ijk = i_j_k_to_ijk(n,I,J,k);
          fprintf(file," %0.15f ",patch->node[ijk]->X[2]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      
      Free(fields);
      Fclose(file);
    }
  }    
  
  free_2d(flds);
  free_2d(d);
  
  FUNC_TOC
}
