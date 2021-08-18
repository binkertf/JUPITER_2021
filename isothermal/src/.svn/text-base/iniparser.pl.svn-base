#!/usr/bin/perl
open INI, "libinit.txt";
@ini = <INI>;
close INI;
open INIC, ">libinit.cx";
open CODE, ">initcode.c";
open DUMP, ">dumpinit.c";
open INC, ">includes/init.h";
print INIC <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "iniparser.pl" for details   */
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
/* "iniparser.pl" for details   */
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
/* "iniparser.pl" for details   */
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
/* "iniparser.pl" for details   */
/*                              */
/********************************/

#include "init.h"
#include "initcodes.h"
#include "jupiter.h"

void DumpInitCode (filename)
char *filename;
{
    FILE *init;
    long i;
    init = prs_open (filename);
    for (i = 0; i < NbFluids; i++) {
    findcodeinit (InitCodeNames[i], &InitCode);
    fprintf (init, "Fluid %ld (%s):\\n", i, FluidName[i]);
    switch (InitCode) {
EOF
print "Running initial conditions parser\n";
$firstini = 1;
$counter = 0;
print INIC "testdoubledefined (InitCode, density);\n";
print INIC "switch (InitCode) {\n";
print CODE "#include \"init.h\"\n#include \"initcodes.h\"\n";
print CODE "#include \"jupiter.h\"\n";
print CODE "static boolean InitCodesSet = NO;\n";
print CODE "void setinitcodes () {\n";
print CODE "if (InitCodesSet) return;\n";
print CODE "InitCodesSet = YES;\n";
print INC "#define undef 0\n";
foreach (@ini) {
  if (($_ !~ /^\s*\#/) && ($_ !~ /^\s$/)) {
    if (/^(\w+)(:|)/) {
      if (!$firstini) {print INIC "\tbreak;\n"; print DUMP "\tbreak;\n"}
      print INIC "case ".$1.":\n";
      print DUMP "case ".$1.":\n";
      $counter++;
      print INC "#define ".$1."\t".$counter."\n";
      print CODE "addinitcode (\"".$1."\", ".$1.");\n";
      $firstini = 0;
    }
    if (/^\s+(.*)$/) {
      print INIC "\t".$1."\n";
      print DUMP "\tfprintf (init, \"   ".$1."\\n\");\n";
    }
  }
}
print INIC "\tbreak;\ndefault:\n\ttestundef (density);\n\tbreak;\n}\n";
print DUMP "\tbreak;\ndefault:\n";
open INIT, "Init.c";
@init = <INIT>;
close INIT;
$exec_write = 0;
foreach (@init) {
    if ($_ =~ 	/between this comment and the next one/) {$exec_write = 1;}
    if ($_ =~ 	/above this comment/) {$exec_write = 0;}
    if (($_ !~ /^\s*\/\*/) &&  ($_ !~ /^\s$/) && ($exec_write)) {
	chomp $_;
	print DUMP "\tfprintf (init, \"".$_."\\n\");\n";
    }
}
print DUMP "\tbreak;\n}\n}\nfclose (init);\n}\n";
print CODE "}\n";
close INIC;
close INC;
close CODE;
