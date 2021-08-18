#include "jupiter.h"
#define POTCODE
#include "potcodes.h"
#include "pot.h"

static long CodeNb;

void addpotcode (string, code)
     char *string;
     long code;
{
  strcpy (PotLib[CodeNb].codename, string);
  PotLib[CodeNb].index = code;
  CodeNb++;
}

void testpotdoubledefined (potcode, pot)
     long potcode;
     real pot;
{
  if ((potcode > 0) && (pot > -1e10)) {
    prs_stderr ("Potential expression defined twice\n");
    prs_stderr ("You should remove any definition in Extpot.c\n");
    prs_error  ("Or you should remove the POTENTIALCODE parameter.");
  }
}

void testundefpot (potential)
     real potential;
{
  if (potential < -1e19) {
    prs_stderr ("You do not have specified the external potential.\n");
    prs_stderr ("Either you edit Extpot.c to define it manually,\n");
    prs_stderr ("or you specify a valid POTENTIALCODE parameter\n");
    prs_error  ("(see src/libpot.txt for a list of these codes)");
  }
}

void findcodepot (string, code)
     char *string;
     long *code;
{
  long i=0;
  setpotcodes ();
  if (*string != 0) {
    i = 0;
    while (strcmp (string, PotLib[i].codename)) {
      if (i >= CodeNb)
	prs_error ("Unknow Potential code");
      i++;
    } 
    *code = i+1;
  } else *code = 0;
}
  
