#include "jupiter.h"

void RT_main (real dt) {
  ComputeTemperatureField ();
  if (NbFluids>1){
    ComputeOpacity ();
  }else{
    ComputeOpacityInit ();
  }
  ComputeStellarHeating();
  ComputeRadiativeEnergy (dt);
  ComputeDiffusionCoefficients ();
  ComputeMatrixElements (dt);
  PredictNewRadEnergy (dt);
  SetRTBoundaryConditions ();
  SolveMatrix (NO); // Solve elliptic equation on non-refined zones only
  GetRadEnergyDifference (dt);
  RadEnergyToTemperature ();
  TemperatureToEnergy ();
}
