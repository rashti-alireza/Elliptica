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
// note: parameter "txt_output_1d" supports regular expression. 
// note: the user must provide the coordinate in the reference interval,
// i.e.,[0,1]x[0,1]x[0,1]. this makes more sense where we are interested 
// to plot the line irrespective of patch coordinates. 
// a linear map then changes this interval according to the patch 
// reference coords (X,Y,Z). */
void print_fields_1D(const Grid_T *const grid,const int iteration, 
                      const char *const folder)
{
  FUNC_TIC
  printf(Pretty0"iteration = %d\n",iteration);
  
  /* pars struct for the print line */
  struct pars_s
  {
    Uint Xline:1;/* if X-line is asked 1, otherwise 0. */
    Uint Yline:1;/* if Y-line is asked 1, otherwise 0. */
    Uint Zline:1;/* if Z-line is asked 1, otherwise 0. */
    double X;/* X coords of the print line */
    double Y;/* Y coords of the print line */
    double Z;/* Z coords of the print line */
    char strv[STR_LEN/2];/* save the line str given in the param file. */
    char suffix[STR_LEN/2];/* suffix for the file. */
  };
  const char *const line_par_name  = "txt_output_1d_line";
  const char *const field_par_name = "txt_output_1d";
  /* the user must provide the coordinate in the reference interval, i.e.,
  // [0,1]x[0,1]x[0,1]. this makes more sense where we are interested 
  // to plot the line irrespective of patch coordinates. 
  // NOTE: if you change this, please change the error msg too. */
  const double REF_coord_min[3] = {0,0,0};
  const double REF_coord_max[3] = {1,1,1};
  const int map_type = 0;/* for future if you want to change the map */
  /* list of the fields(flds) to be printed out */
  char **flds = 
     read_separated_items_in_string(PgetsEZ(field_par_name),',');
  /* list of lines(lns) */
  char **lns = 
     read_separated_items_in_string(PgetsEZ(line_par_name),')');
  Uint Nlines;/* number of lines (pars struct) */
  Uint l,p,counter;
  
  if (!flds)
    printf(Pretty0"No field is specified.\n");
  
  if (!lns)
    printf(Pretty0"No line is specified.\n");
  
  /* count the number of lines to create the struct */
  counter = 0;
  for (l = 0; lns && lns[l]; ++l)
  {
    if (strlen(lns[l]) > 5)/* coz we need, e.g., "(I,1,2" and not ",".
                           // note ')' is eliminated. */
      counter++;
  }
  Nlines = counter;
  struct pars_s parsed[counter];
  
  /* populate the parsed structs */
  p = 0;/* index for parsed */
  for (l = 0; lns && lns[l]; ++l)
  {
     if (strlen(lns[l]) <= 5) continue;
     
     char **subs = 0;
     
     /* substitue ,e.g., ",(0,J,2" by "0,J,2". */
     regex_replace(lns[l],"^.?\\(","",lns[l]);
     subs = read_separated_items_in_string(lns[l],',');
     
     /* sanity checks: */
     /* nothing found */
     if (!subs) Errors("Wrong format for %s.\n",line_par_name);
     /* we need 3 pieces like 0,J, and 2  in "(0,J,2)" */
     for (counter = 0; subs[counter]; counter++);
     if (counter != 3) Errors("Wrong format for %s.\n",line_par_name);
     /* correct position */
     if (
         /* 'I/i' should only take place on the 0th position */
         *subs[1] == 'X' || 
         *subs[1] == 'x' || 
         *subs[2] == 'X' || 
         *subs[2] == 'x' ||
         /* 'J/j' should only take place on the 1st position */
         *subs[0] == 'Y' || 
         *subs[0] == 'y' || 
         *subs[2] == 'Y' || 
         *subs[2] == 'y' || 
         /* 'K/k' should only take place on the 2nd position */
         *subs[0] == 'Z' || 
         *subs[0] == 'z' || 
         *subs[1] == 'Z' || 
         *subs[1] == 'z'
         )
     {
       Errors("Wrong order for %s.\n",line_par_name); 
     }
     
     /* cp str */
     sprintf(parsed[p].strv,"(%s)",lns[l]);
     
     /* realize the components */
     counter = 0;
     if (*subs[0] == 'X' || *subs[0] == 'x')
     {
       parsed[p].Xline = 1;
       parsed[p].Yline = 0;
       parsed[p].Zline = 0;
       
       parsed[p].X = DBL_MAX;
       parsed[p].Y = atof(subs[1]);
       parsed[p].Z = atof(subs[2]);
       sprintf(parsed[p].suffix,"X_%s_%s",subs[1],subs[2]);
       
       /* interval check */
       if (LSS(parsed[p].Y,REF_coord_min[1]) ||
           GRT(parsed[p].Y,REF_coord_max[1]) ||
           LSS(parsed[p].Z,REF_coord_min[2]) ||
           GRT(parsed[p].Z,REF_coord_max[2]) 
           )
           Errors("%s falling outside of the reference interval [0,1]x[0,1]x[0,1].\n",
                  parsed[p].strv);
       counter++;
     }
     if (*subs[1] == 'Y' || *subs[1] == 'y')
     {
       parsed[p].Xline = 0;
       parsed[p].Yline = 1;
       parsed[p].Zline = 0;
       
       parsed[p].X = atof(subs[0]);
       parsed[p].Y = DBL_MAX;
       parsed[p].Z = atof(subs[2]);
       sprintf(parsed[p].suffix,"%s_Y_%s",subs[0],subs[2]);
       /* interval check */
       if (LSS(parsed[p].X,REF_coord_min[0]) ||
           GRT(parsed[p].X,REF_coord_max[0]) ||
           LSS(parsed[p].Z,REF_coord_min[2]) ||
           GRT(parsed[p].Z,REF_coord_max[2]) 
           )
           Errors("%s falling outside of the reference interval [0,1]x[0,1]x[0,1].\n",
                  parsed[p].strv);
       counter++;
     }
     if (*subs[2] == 'Z' || *subs[2] == 'z')
     {
       parsed[p].Xline = 0;
       parsed[p].Yline = 0;
       parsed[p].Zline = 1;
       
       parsed[p].X = atof(subs[0]);
       parsed[p].Y = atof(subs[1]);
       parsed[p].Z = DBL_MAX;
       sprintf(parsed[p].suffix,"%s_%s_Z",subs[0],subs[1]);
       /* interval check */
       if (LSS(parsed[p].Y,REF_coord_min[1]) ||
           GRT(parsed[p].Y,REF_coord_max[1]) ||
           LSS(parsed[p].X,REF_coord_min[0]) ||
           GRT(parsed[p].X,REF_coord_max[0]) 
           )
           Errors("%s falling outside of the reference interval [0,1]x[0,1]x[0,1].\n",
                  parsed[p].strv);
       counter++;
     }
     /* if we have more than two letters, e.g., "(I,J,7)". */
     if (counter != 1)
       Errors("Wrong format for %s.\n",line_par_name); 
       
     p++;
     free_2d(subs);
  }
  
  /* test */
  if (0)
  for (l = 0; l < Nlines; ++l)
  {
    printf("param[%u] = %s => (%u,%u,%u) & (%g,%g,%g)| %s\n",
            l, parsed[l].strv, 
            parsed[l].Xline,parsed[l].Yline,parsed[l].Zline,
            parsed[l].X,parsed[l].Y,parsed[l].Z,parsed[l].suffix);
  }
  /* check for redundancy */
  for (l = 0; l < Nlines; ++l)
  {
    struct pars_s *line = &parsed[l];
    for (Uint lp = l+1; lp < Nlines; ++lp)
    {
      struct pars_s *linep = &parsed[lp];
      if (!strcmp(line->strv,linep->strv))
        Errors("Multiple line with same value for %s.",line->strv);
    }
  }

  if(flds && lns)
  for (l = 0; l < Nlines; ++l)
  {
    struct pars_s *line = &parsed[l];
    
    FOR_ALL_PATCHES(p,grid)
    {
      Patch_T *patch = grid->patch[p];
      const Uint *const n = patch->n;
      const char *stem = strstr(patch->name,"_"); stem++;
      Field_T **fields = 0;/* save the found fields */
      Uint Nfld      = 0;/* total num of the found fields. */
      FILE *file     = 0;
      char file_name[STR_LEN] = {0};
      double min_d;
      double X,Y,Z;
      Uint I,J,K;
      Uint i,j,k,ijk,f;
      
      /* finding the index of the closest point to 
      // the specified one in the param. */
      I = J = K = UINT_MAX;/* to catch bug */
      if (line->Xline)
      {
        /* normalize */
        Y = map_to_ref_interval(line->Y,patch,REF_coord_min,REF_coord_max,1,map_type);
        Z = map_to_ref_interval(line->Z,patch,REF_coord_min,REF_coord_max,2,map_type);
       
        min_d = DBL_MAX;
        for (j = 0; j < n[1]; ++j)
        {
          ijk = i_j_k_to_ijk(n,0,j,0);
          double dist = fabs(Y-patch->node[ijk]->X[1]);
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
          double dist = fabs(Z-patch->node[ijk]->X[2]);
          if (dist < min_d)
          {
            K = k;
            min_d = dist;
          }
        }
      }
      
      else if (line->Yline)
      {
        /* normalize */
        X = map_to_ref_interval(line->X,patch,REF_coord_min,REF_coord_max,0,map_type);
        Z = map_to_ref_interval(line->Z,patch,REF_coord_min,REF_coord_max,2,map_type);

        min_d = DBL_MAX;
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,0,0);
          double dist = fabs(X-patch->node[ijk]->X[0]);
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
          double dist = fabs(Z-patch->node[ijk]->X[2]);
          if (dist < min_d)
          {
            K = k;
            min_d = dist;
          }
        }
      }
      
      else if (line->Zline)
      {
        /* normalize */
        X = map_to_ref_interval(line->X,patch,REF_coord_min,REF_coord_max,0,map_type);
        Y = map_to_ref_interval(line->Y,patch,REF_coord_min,REF_coord_max,1,map_type);

        min_d = DBL_MAX;
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,0,0);
          double dist = fabs(X-patch->node[ijk]->X[0]);
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
          double dist = fabs(Y-patch->node[ijk]->X[1]);
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
      fields = find_field_by_name_or_regex(patch,flds,&Nfld);
      
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
        fprintf(file,"# line_coordinate x(X,Y,Z) y(X,Y,Z) z(X,Y,Z)");
        for (f = 0; f < Nfld; ++f)
          fprintf(file," %s", fields[f]->name);
        fprintf(file,"\n");
      }
      fprintf(file,"\n# \"time = %d\", (X,Y,Z) = %s,",iteration,line->strv);
      
      /* print coords and fields in each column. */
      if (line->Xline)
      {
        Y = patch->node[i_j_k_to_ijk(n,0,J,0)]->X[1];
        Z = patch->node[i_j_k_to_ijk(n,0,0,K)]->X[2];
        fprintf(file," Y[%u] = %g, Z[%u] = %g\n",J,Y,K,Z);
        
        for (i = 0; i < n[0]; ++i)
        {
          ijk = i_j_k_to_ijk(n,i,J,K);
          fprintf(file,"%0.15f %0.15f %0.15f %0.15f",
                        patch->node[ijk]->X[0],
                        patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      else if (line->Yline)
      {
        X = patch->node[i_j_k_to_ijk(n,I,0,0)]->X[0];
        Z = patch->node[i_j_k_to_ijk(n,0,0,K)]->X[2];
        fprintf(file," X[%u] = %g, Z[%u] = %g\n",I,X,K,Z);
        
        for (j = 0; j < n[1]; ++j)
        {
          ijk = i_j_k_to_ijk(n,I,j,K);
          fprintf(file,"%0.15f %0.15f %0.15f %0.15f",
                        patch->node[ijk]->X[1],
                        patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      else if (line->Zline)
      {
        X = patch->node[i_j_k_to_ijk(n,I,0,0)]->X[0];
        Y = patch->node[i_j_k_to_ijk(n,0,J,0)]->X[1];
        fprintf(file," X[%u] = %g, Y[%u] = %g\n",I,X,J,Y);
        
        for (k = 0; k < n[2]; ++k)
        {
          ijk = i_j_k_to_ijk(n,I,J,k);
          fprintf(file,"%0.15f %0.15f %0.15f %0.15f",
                        patch->node[ijk]->X[2],
                        patch->node[ijk]->x[0],
                        patch->node[ijk]->x[1],
                        patch->node[ijk]->x[2]);
          for (f = 0; f < Nfld; ++f)
            fprintf(file," %0.15f",fields[f]->v[ijk]);
          fprintf(file,"\n");
        }
      }
      
      Free(fields);
      Fclose(file);
    }/* end of FOR_ALL_PATCHES(p,grid) */
  }    
  
  free_2d(flds);
  free_2d(lns);
  
  FUNC_TOC
}

/* :-> a set of pointers to the found fields in the given patch.
// the number of pointers is set in Nfld and the field names 
// given should be given in the argument fld_names. */
static Field_T **find_field_by_name_or_regex(const Patch_T *const patch,
                                             char **const fld_names,
                                             Uint *const Nfld)
{
  Field_T **fields = 0;
  *Nfld = 0;
  
  for (Uint f = 0; fld_names[f]; ++f)
  {
    int field_ind  = 0;
    Uint Nf        = 0;
    Uint *fInd     = 0;
    Uint Ufield_ind= 0;
    Field_T *field = 0;
    Uint fi;
  
    field_ind  = _Ind(fld_names[f]); 
    /* if couldn't find, maybe its a regex. */
    if (field_ind < 0)
    {
      fInd = find_field_index_regex(patch,fld_names[f],&Nf);
        
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
      
      fields = realloc(fields,(*Nfld+1)*sizeof(*fields));
      IsNull(fields);
      fields[*Nfld] = field;
      (*Nfld)++;
    }
    /* only if it is a regex found free it. */
    if (field_ind < 0) Free(fInd);
  }/* end of for (f = 0; fld_names[f]; ++f) */
  
  return fields;
}

/* :-> map the given point x to the [patch->min,patch->max] interval. */
static double map_to_ref_interval(const double X,
                                  const Patch_T *const patch,
                                  const double *const min,
                                  const double *const max,
                                  const int dir, const int map_type)
{
  const double mp = min[dir];
  const double Mp = max[dir];
  const double m  = patch->min[dir];
  const double M  = patch->max[dir];
  double a,b;
  double y = DBL_MAX;
  
  switch (map_type)
  {
    case 0:/* linear map f: X-> aX+b */
      
      a = (m - M)/(mp - Mp);
      b = (M*mp - m*Mp)/(mp - Mp);
      y = a*X+b;
      
      break;
    default:
      Error0(NO_OPTION);
  }
  
  /* test */
  //printf("dir[%d]: %g -> %g, [%g,%g]\n",dir,X,y,m,M);
  
  return y;
}
