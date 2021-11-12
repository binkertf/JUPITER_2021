#define __VAR_DEF
#include "jupiter.h"
#undef __VAR_DEF

void InitVariables ()
{
  var("DT", &DT, REAL, YES, "1.");
  var("OUTPUTDIR", OUTPUTDIR, STRING, YES, "");
  var("NTOT", &NTOT, INT, YES, "64.0");
  var("NINTERM", &NINTERM, INT, YES, "64.0");
  var("SIZE1", &SIZE1, INT, YES, "64.0");
  var("SIZE2", &SIZE2, INT, YES, "64.0");
  var("SIZE3", &SIZE3, INT, YES, "64.0");
  var("GHOSTFILLINGORDER", &GHOSTFILLINGORDER, INT, NO, "1.0");
  var("DIM1SPACING", DIM1SPACING, STRING, NO, "Arithmetic");
  var("DIM2SPACING", DIM2SPACING, STRING, NO, "Arithmetic");
  var("DIM3SPACING", DIM3SPACING, STRING, NO, "Arithmetic");
  var("RANGE1LOW", &RANGE1LOW, REAL, NO, "0.0");
  var("RANGE1HIGH", &RANGE1HIGH, REAL, NO, "1.0");
  var("RANGE2LOW", &RANGE2LOW, REAL, NO, "0.0");
  var("RANGE2HIGH", &RANGE2HIGH, REAL, NO, "1.0");
  var("RANGE3LOW", &RANGE3LOW, REAL, NO, "0.0");
  var("RANGE3HIGH", &RANGE3HIGH, REAL, NO, "1.0");
  var("COURANTNUMBER", &CFLSECURITY, REAL, NO, "0.9");
  var("GAMMA", &GAMMA, REAL, NO, "1.666666666");
  var("NDIM", &NDIM, INT, NO, "3");
  var("COORDTYPE", COORDTYPE, STRING, NO, "CARTESIAN");
  var("COORDPERMUT", COORDPERMUT, STRING, NO, "123");
  var("INHIBREFDIM1", &INHIBREFDIM1, BOOL, NO, "NO");
  var("INHIBREFDIM2", &INHIBREFDIM2, BOOL, NO, "NO");
  var("INHIBREFDIM3", &INHIBREFDIM3, BOOL, NO, "NO");
  var("DIM1PERIODIC", &DIM1PERIODIC, BOOL, NO, "NO");
  var("DIM2PERIODIC", &DIM2PERIODIC, BOOL, NO, "NO");
  var("DIM3PERIODIC", &DIM3PERIODIC, BOOL, NO, "NO");
  var("GRIDFILE", GRIDFILE, STRING, NO, "");
  var("RIEMANNSOLVER", RIEMANNSOLVER, STRING, NO, "2R");
  var("DUSTSOLVER", DUSTSOLVER, STRING, NO, "FO");
  var("ADIABATIC", ADIABATIC, STRING, NO, "NO");
  var("METHOD", METHOD, STRING, NO, "MUSCL");
  var("INITCODE", INITCODE, STRING, NO, "");
  var("POTENTIALCODE", POTENTIALCODE, STRING, NO, "");
  var("CONDLIM1", CONDLIM1, STRING, NO, "UNDEFINED");
  var("CONDLIM2", CONDLIM2, STRING, NO, "UNDEFINED");
  var("CONDLIM3", CONDLIM3, STRING, NO, "UNDEFINED");
  var("FIXEDPOTENTIAL", FIXEDPOTENTIAL, STRING, NO, "NEVER");
  var("KEPLERIAN", &KEPLERIAN, BOOL, NO, "NO");
  var("OMEGAFRAME", &OMEGAFRAME, REAL, NO, "0.0");
  var("ASPECTRATIO", &ASPECTRATIO, REAL, NO, "0.05");
  var("SIGMASLOPE", &SIGMASLOPE, REAL, NO, "0.0");
  var("FLARINGINDEX", &FLARINGINDEX, REAL, NO, "0.0");
  var("SIGMA0", &SIGMA0, REAL, NO, "6e-4");
  var("MASSTAPER", &MASSTAPER, REAL, NO, "0.0");
  var("VISCOSITY", &VISCOSITY, REAL, NO, "0.0");
  var("PLANETMASS", &PLANETMASS, REAL, NO, "0.0");
  var("SMOOTHING", &SMOOTHING, REAL, NO, "0.0");
  var("GRIDFRICTION1", &GRIDFRICTION1, REAL, NO, "0.0");
  var("GRIDFRICTION2", &GRIDFRICTION2, REAL, NO, "0.0");
  var("GRIDFRICTION3", &GRIDFRICTION3, REAL, NO, "0.0");
  var("EXTERNALPOTENTIAL", &EXTERNALPOTENTIAL, BOOL, NO, "NO");
  var("MERIDIANCORRECTION", &MERIDIANCORRECTION, BOOL, NO, "NO");
  var("DRIFTVELOCITY", &DRIFTVELOCITY, REAL, NO, "0.0");
  var("CS", &CS, REAL, NO, "1.0");
  var("INCLINATION", &INCLINATION, REAL, NO, "0.0");
  var("SUBCYCLING", SUBCYCLING, STRING, NO, "Auto");
  var("NODAMPING", NODAMPING, STRING, NO, "NO");
  var("POTENTIALHYDROSTAT", POTENTIALHYDROSTAT, STRING, NO, "");
  var("INITHYDROSTAT", INITHYDROSTAT, STRING, NO, "");
  var("SYMDIMHYDROSTAT", SYMDIM_EQ, STRING, NO, "000");
  var("CORRDIMHYDROSTAT", CORRDIM_EQ, STRING, NO, "000");
  var("FLUIDS", FLUIDS, STRING, NO, "gas");
  var("COUPLING", &COUPLING, REAL, NO, "0.0");
  var("DUSTSIZE", &DUSTSIZE, REAL, NO, "0.001");
  var("DUSTSOLIDRHO", &DUSTSOLIDRHO, REAL, NO, "3.0");
  var("DUSTTOGAS", &DUSTTOGAS, REAL, NO, "0.01");
  var("DUSTDENSFLOOR", &DUSTDENSFLOOR, REAL, NO, "1e-15");
  var("OUTPUTATREFINEMENT", &OUTPUTATREFINEMENT, INT, NO, "0");
  var("FINESTLEVEL", &FinestLevel, INT, NO, "0");
  var("HIGHRESLEVEL", &HIGHRESLEVEL, INT, NO, "100");
  var("VISCUTOFFLEVEL", &VISCUTOFFLEVEL, INT, NO, "100");
  var("ENFORCEETOTCONSERVATION", &ETOTCONS, BOOL, NO, "YES");
  var("PRESSURERATIO", &PRESSURERATIO, REAL, NO, "1.0");
  var("DISKSPEEDFACTOR", &DISKSPEEDFACTOR, REAL, NO, "1.0");
  var("HALFDISK", HALFDISK, STRING, NO, "YES");
  var("RADIATIVE", RADIATIVE, STRING, NO, "NO");
  var("STELLAR", STELLAR, STRING, NO, "NO");
  var("TSTAR", &TSTAR, REAL, NO, "5600.0");
  var("RSTAR", &RSTAR, REAL, NO, "3.0");
  var("STELLDEG", &STELLDEG, REAL, NO, "7.0");
  var("SMRATIO", &SMRATIO, REAL, NO, "2.0");
  var("R0", &R0, REAL, YES, "77792000000000.0");
  var("XMSTAR", &XMSTAR, REAL, YES, "1.0");
}
