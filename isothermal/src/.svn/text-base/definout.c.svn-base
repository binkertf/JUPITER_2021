#include "jupiter.h"

/* In this file we read the default root directories for output and
 input. They should be defined by the environment variables
 JUPITER_OUT and JUPITER_IN. If these variables are undefined, they
 are set by default to ./out and ./in */

void ReadDefaultInOut () {
  char *defin;
  char *defout;
  defin = getenv ("JUPITER_IN");
  if (defin != NULL) {
    strncpy (DefaultIn, defin, MAXLINELENGTH-2);
  } else {
    sprintf (DefaultIn, "./in");
  }
  if (*(DefaultIn+strlen(DefaultIn)-1) != '/')
    strcat (DefaultIn, "/");	/* Add trailing slash if missing */
  pInfo ("The default input directory root is %s\n", DefaultIn);
  defout = getenv ("JUPITER_OUT");
  if (defout != NULL) {
    strncpy (DefaultOut, defout, MAXLINELENGTH-2);
  } else {
    sprintf (DefaultOut, "./out");
  }
  if (*(DefaultOut+strlen(DefaultOut)-1) != '/')
    strcat (DefaultOut, "/");	/* Add trailing slash if missing */
  pInfo ("The default output directory root is %s\n", DefaultOut);
}

void SubsDef (msg, target, def)
  char *msg, *target, *def;
{
  char c='@';
  char new_target[MAXLINELENGTH];
  char *loc, *follow;
  loc = strchr (target, (int)c);
  if (loc != NULL) {
    if (*(loc+1) == '/')
      follow=loc+2;
    else
      follow=loc+1;
    if (loc != target) {
      pWarning ("Characters located before '%c' wildcard in %s definition are discarded\n", c, msg);
    }
    snprintf (new_target,MAXLINELENGTH/2-1,def);
    strncat (new_target,follow,MAXLINELENGTH/2-1);
    strncpy (target, new_target, MAXLINELENGTH-1);
  }
  pInfo ("The %s is %s\n", msg, target);
  prs_msg ("The %s is %s\n", msg, target);
}

void ExtractPath (filename, path, shortname)
     char *filename, *path, *shortname;
{
  char *lastslash;
  /* We assume that filename and path have MAXLINELENGTH allocated chars */
  strncpy (path, filename, MAXLINELENGTH-1);
  lastslash = strrchr (path, '/');
  if (lastslash == NULL) {
    sprintf (path, "./");
    strcpy (shortname, filename);
  }
  else {
    strcpy (shortname, lastslash+1);
    *(lastslash+1) = 0;
  }
}
