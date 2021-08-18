#!/usr/bin/perl
open INI, "libpot.txt";
@ini = <INI>;
close INI;
open INIC, ">libpot.cx";
open CODE, ">potcode.c";
open DUMP, ">dumppot.c";
open INC, ">includes/pot.h";
print INIC <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "potparser.pl" for details   */
/*                              */
/********************************/
EOF
print INC <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "potparser.pl" for details   */
/*                              */
/********************************/
EOF
print CODE <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "potparser.pl" for details   */
/*                              */
/********************************/
EOF
print DUMP <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "potparser.pl" for details   */
/*                              */
/********************************/

#include "pot.h"
#include "potcodes.h"
#include "jupiter.h"

void DumpPotCode (filename)
char *filename;
{
    FILE *pot;
    pot = prs_open (filename);
    if (!EXTERNALPOTENTIAL) {
	fprintf (pot, "No external potential is applied to the system.\\n");
	fclose (pot);
	return;
    }
    switch (PotentialCode) {
EOF
print "Running external potential parser\n";
$firstini = 1;
$counter = 0;
print INIC "if (NeverTested) {\n";
print INIC "    testpotdoubledefined (PotentialCode, pot);\n";
print INIC "    NeverTested = FALSE;\n";
print INIC "}\n";
print INIC "switch (PotentialCode) {\n";
print CODE "#include \"pot.h\"\n#include \"potcodes.h\"\n";
print CODE "#include \"jupiter.h\"\n";
print CODE "static boolean PotCodesSet = NO;\n";
print CODE "void setpotcodes () {\n";
print CODE "if (PotCodesSet) return;\n";
print CODE "PotCodesSet = YES;\n";
foreach (@ini) {
  if (($_ !~ /^\s*\#/) && ($_ !~ /^\s$/)) {
    if (/^(\w+)(:|)/) {
      if (!$firstini) {print INIC "\tbreak;\n"; print DUMP "\tbreak;\n";}
      print INIC "case ".$1.":\n";
      print DUMP "case ".$1.":\n";
      $counter++;
      print INC "#define ".$1."\t".$counter."\n";
      print CODE "addpotcode (\"".$1."\", ".$1.");\n";
      $firstini = 0;
    }
    if (/^\s+(.*)$/) {
      print INIC "\t".$1."\n";
      print DUMP "\tfprintf (pot, \"".$1."\\n\");\n";
    }
  }
}
print DUMP "\tbreak;\ndefault:\n";
open POT, "ExtPot.c";
@pot = <POT>;
close POT;
$exec_write = 0;
foreach (@pot) {
    if ($_ =~ 	/remove the present potential set up/) {$exec_write = 1;}
    if ($_ =~ 	/Use any of the following/) {$exec_write = 0;}
    chomp $_;
    if (($_ !~ /^\s*\/\*/) &&  ($_ !~ /^\s*$/) && ($exec_write)) {
	print DUMP "\tfprintf (pot, \"".$_."\\n\");\n";
    }
}
print DUMP "\tbreak;\n}\nfclose (pot);\n}\n";
print INIC "\tbreak;\ndefault:\n\ttestundefpot (pot);\n\tbreak;\n}\n";
print CODE "}\n";
close INIC;
close INC;
close CODE;
