#include "jupiter.h"

extern char *SCoordNames[];

void GridFileError (long linenumber, char *filename, const char *template, ...)
{
  char string[MAXLINELENGTH];
  va_list ap;
  va_start (ap, template);
  vsprintf (string, template, ap);
  va_end (ap);
  prs_error ("%s in %s at line %ld.\n", string, filename, linenumber);
}

FILE *FindGridFile (filename)
     char *filename;
{
  char cwd[MAXLINELENGTH];
  char loc1[MAXLINELENGTH];
  char loc2[MAXLINELENGTH];
  char loc3[MAXLINELENGTH];
  char finalname[MAXLINELENGTH];
  FILE *grid;
  pInfo ("The grid file name is %s before expansion\n", filename);
  SubsDef ("grid file", filename, DefaultIn);
  if (strchr (filename, '/') == NULL) {
    pInfo ("No path is specified.\n");
    pInfo ("Therefore I look for this file at the following locations:\n");
    pInfo ("1/ Directory where the parameter file resides (%s)\n", InputPath);
    pInfo ("2/ Current working directory (%s)\n", getcwd(cwd,MAXLINELENGTH-1));
    pInfo ("Probing location 1/\n");
    sprintf (finalname, "%s%s", InputPath, filename);
    pInfo ("Does %s exist ? ", finalname);
    grid = fopen (finalname, "r");
    if (grid != NULL) {
      pInfo ("YES\nUsing grid file %s\n", finalname);
      return grid;
    } else {
      pInfo ("NO\nProbing location 2/\n");
      sprintf (finalname, "%s/%s", cwd, filename);
      pInfo ("Does %s exist ? ", finalname);
      grid = fopen (finalname, "r");
      if (grid != NULL) {
	pInfo ("YES\nUsing grid file %s\n", finalname);
	return grid;
      } else {
	pInfo ("NO\nGrid file not found neither at location 1 nor 2, giving up.\n");
	prs_error ("Grid file %s not found in  %s nor in %s. Please check location.\n", \
		filename, InputPath, cwd);
      }
    }
  } else {
    /* Filename contains a slash */
    if (*filename == '/') {
      pInfo ("An absolute path is specified.\n");
      grid = fopen (filename, "r");
      if (grid != NULL) {
	pInfo ("The grid file %s exists. OK\n", filename);
	return grid;
      } else {
	pInfo ("The grid file %s does not exist. Aborting.\n", filename);
	prs_error ("The grid file %s does not exist. Please check location.\n", filename);
      }
    } else {
      pInfo ("A relative path is specified.\n");
      pInfo ("Therefore I look for this file at the following locations:\n");
      pInfo ("1/ Subdirectory of current working directory (%s)\n", getcwd(cwd,MAXLINELENGTH-1));
      pInfo ("2/ Subdirectory of the directory where the parameter file resides (%s)\n", InputPath);
      pInfo ("3/ Subdirectory of the root input file directory (%s), in case you forgot an '@' wildcard.\n", DefaultIn);
      pInfo ("Probing location 1/\n");
      sprintf (finalname, "%s/%s", cwd, filename);
      pInfo ("Does %s exist ? ", finalname);
      strcpy (loc1, finalname);
      grid = fopen (finalname, "r");
      if (grid != NULL) {
	pInfo ("YES\nUsing grid file %s\n", finalname);
	return grid;
      } else {
	pInfo ("NO\nProbing location 2/\n");
	sprintf (finalname, "%s%s", InputPath, filename);
	pInfo ("Does %s exist ? ", finalname);
	strcpy (loc2, finalname);
	grid = fopen (finalname, "r");
	if (grid != NULL) {
	  pInfo ("YES\nUsing grid file %s\n", finalname);
	  return grid;
	} else {
	  pInfo ("NO\nProbing location 3/\n");
	  sprintf (finalname, "%s%s", DefaultIn, filename);
	  pInfo ("Does %s exist ? ", finalname);
	  strcpy (loc3, finalname);
	  grid = fopen (finalname, "r");
	  if (grid != NULL) {
	    pInfo ("YES\nUsing grid file %s\n", finalname);
	    return grid;
	  } else {
	    pInfo ("NO\nGrid file not found at locations 1, 2 and 3. Giving up.\n");
	    prs_error ("I could not locate the grid file in\n%s,\n%s, or\n%s\nPlease check location.\n", loc1, loc2, loc3);
	  }
	}
      }
    }
  }
  return NULL;
}

void ScanGridFile (filename)
     char *filename;
{
  FILE *input, *output;
  char *shortname;
  GridFileInfo grids[MAXGRIDS];
  long line=0, bc[6]={0,0,0,0,0,0}, nb = 0;
  double xmin[3], xmax[3], foo;	/* 'foo' is used to check trailing parameters */
  char string[MAXLINELENGTH];
  static long number=0, expect, level, i, j;
  long linenumber=0;
  input = FindGridFile (filename);
  if (!(shortname = strrchr (filename, '/'))) shortname = filename-1;
  if (!CPU_Rank) {
    output = prs_open (++shortname);
    switch (NDIM) {
    case 1: 
      fprintf (output, "# %s_min %s_max level BC_%s_min BC_%s_max\n\n",\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]]);
      break;
    case 2: 
      fprintf (output, "# %s_min %s_min %s_max %s_max level BC_%s_min BC_%s_min BC_%s_max BC_%s_max\n\n",\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]]);
      break;
    case 3: 
      fprintf (output, "# %s_min %s_min %s_min %s_max %s_max %s_max level BC_%s_min BC_%s_min BC_%s_min BC_%s_max BC_%s_max BC_%s_max\n\n",\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[2]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[2]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[2]],\
	       SCoordNames[CoordType*3+InvCoordNb[0]],\
	       SCoordNames[CoordType*3+InvCoordNb[1]],\
	       SCoordNames[CoordType*3+InvCoordNb[2]]);
      break;
    }
  }
  expect = 1+4*NDIM;
  while (!feof(input)) {
    if (fgets (string, MAXLINELENGTH, input) != NULL) {
      linenumber++;
      if ((*string != '#') && (strlen(string) > 1)) {
	switch (NDIM) {
	case 3: 
	  nb = sscanf (string, "%lf %lf %lf %lf %lf %lf %ld %ld %ld %ld %ld %ld %ld %lf", \
		       xmin, xmin+1, xmin+2, xmax, xmax+1, xmax+2, &level,\
		       bc, bc+1, bc+2, bc+3, bc+4, bc+5, &foo);
	  break;
	case 2: 
	  nb = sscanf (string, "%lf %lf %lf %lf %ld %ld %ld %ld %ld %lf", \
		       xmin, xmin+1, xmax, xmax+1, &level,\
		       bc, bc+1, bc+3, bc+4, &foo);
	  break;
	case 1: 
	  nb = sscanf (string, "%lf %lf %ld %ld %ld %lf", \
		       xmin, xmax, &level,\
		       bc, bc+3, &foo);
	  break;
	}
	if (nb < expect)
	  GridFileError (linenumber, filename,\
			 "Badly formed refinement instruction (missing values ?)");
	if (nb > expect)
	  GridFileError (linenumber, filename,\
			 "Badly formed refinement instruction (more than %ld values)",\
			 expect);
	if (level > MaxLevel) {
	  pInfo ("Discarding line %d in gridfile since its level (%d) is higher than %d\n",\
		 linenumber, level, MaxLevel);
	} else {
	  if (!CPU_Rank) fputs (string,  output);
	  grids[line].number = number++;
	  grids[line].linenumber = linenumber;
	  for (i = 0; i < 3; i++)	{ /* 3, not NDIM */
	    grids[line].xmin[i] = (i < NDIM ? xmin[i]: corner_min0[i]);
	    grids[line].xmax[i] = (i < NDIM ? xmax[i]: corner_min0[i]); 
	    /* min as well for trailing dimensions */
	    for (j = INF; j <= SUP; j++) {
	      grids[line].bc[i+j*3] = bc[i+j*3];
	    }
	  }
	  grids[line].level = level;
	  grids[line].last = TRUE;
	  if (line > 0) grids[line-1].last = FALSE;
	  line++;
	}
	if (!(line < MAXGRIDS)) {
	  prs_stderr ("Too many grid specified. You need to recompile the code\n");
	  prs_error  ("after increasing the value of MAXGRIDS in includes/def.h");
	}
      }
    }
  }
  fclose (input);
  if (!CPU_Rank) fclose (output);
  GridAbs (grids);
  GridPos (grids);
  GridBuild (grids);
}
