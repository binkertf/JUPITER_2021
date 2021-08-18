#include "jupiter.h"
#define INITCODE
#include "initcodes.h"
#include "init.h"

static long CodeNb;

void addinitcode (string, code)
     char *string;
     long code;
{
  strcpy (InitLib[CodeNb].codename, string);
  InitLib[CodeNb].index = code;
  CodeNb++;
}

void testdoubledefined (initcode, density)
     long initcode;
     real density;
{
  if ((initcode > 0) && (density > -1e10)) {
    prs_stderr ("Initial conditions defined twice\n");
    prs_stderr ("You should remove any definition in Init.c\n");
    prs_error  ("Or you should remove the INITCODE parameter.");
  }
}

void testundef (density)
     real density;
{
  if (density < 0.0) {
    prs_stderr ("You do not have specified initial conditions.\n");
    prs_stderr ("Either you edit Init.c to define manually IC,\n");
    prs_error  ("or you specify a valid INITCODE parameter.");
  }
}

void findcodeinit (string, code)
     char *string;
     long *code;
{
  long i=0;
  setinitcodes ();
  if (*string != 0) {
    i = 0;
    while (strcmp (string, InitLib[i].codename)) {
      if (i >= CodeNb)
	prs_error ("Unknow Init code");
      i++;
    } 
    *code = i+1;
  } else *code = 0;
}
  
