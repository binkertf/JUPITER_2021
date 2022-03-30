void ComputeTemperatureField (void);
void ComputeOpacity (void);
void ComputeRadiativeEnergy (real);
void ComputeDiffusionCoefficients (void);
void ComputeMatrixElements (real);
void SetRTBoundaryConditions (void);
void SolveMatrix (boolean);
void RadEnergyToTemperature (void);
void TemperatureToEnergy (void);
void prs_error (const char *template, ...);
void prs_end (const char *template, ...);
void prs_stderr (const char *template, ...);
void prs_msg (const char *template, ...);
void *prs_malloc ();
void InitVariables ();
void var();
void ReadVariables ();
void ReadVarFile ();
void ListVariables ();
FILE *prs_open ();
FILE *prs_openi ();
FILE *prs_opend ();
FILE *prs_opena ();
long switches ();
void PrintUsage ();
void MonitorEnergy ();
void FillEnergyTot();
void DumpParameters ();
void DumpSources ();
void ListVariablesIDL ();
void SetGlobalVariables ();
void FillCPUGrid ();
void TestEndian ();
void splitgrid();
void repartition();
void ScanGridFile ();
void primefactors ();
boolean buildprime ();
void GridFileError (long linenumber, char *filename, const char *template, ...);
FluidPatch *CreateFluidPatch ();
ScalarField *CreateScalarField ();
VectorField *CreateVectorField ();
void FreeScalarField ();
void FreeVectorField ();
InterfaceFlux *CreateInterfaceFlux ();
void FreeInterfaceFlux ();
void FreeFluidPatch ();
void WriteWorkArray ();
void WriteVectorField ();
void WriteScalarField ();
void WriteFluid ();
void InitialCondition ();
void ReallocateSpaceForCurrentFluidPatch ();
void ReallocateSpaceForSecondaryFluidPatch ();
void CreateCurrentFluidPatch ();
void CreateSecondaryFluidPatch ();
void SendToCurrent ();
void SendToSecondary ();
void ComputeExternalPotential ();
void ResetFaceFluxesLevel ();
void ResetFaceFluxes ();
void CurrentToPatch ();
real Mass ();
real TotalMass ();
void ApplyBoundaryConditions ();
void ItereLevel ();
real RecursiveIteration ();
void InitWholeHierarchy ();
real Momentum ();
real TotalMomentum ();
void ConservativeUpdate ();
void ConservativeDustUpdate ();
void DiffusionUpdate ();
real ItereRS ();
void GetStar_TWOSHOCKS ();
void GetStar_TWORAREFACTIONS ();
void GetStar_ITERATIVE ();
void GetStar_PRESSURELESS ();
boolean GetStar_ADI_SS ();
boolean GetStar_ADI_RS ();
boolean GetStar_ADI_SR ();
boolean GetStar_AdiabaticSolver ();
real TVDslope();
real minmod();
real superbee();
void FillBeam ();
void FillBeam2 ();
void FillFluxes  ();
void FillDiffFluxes ();
void FillDust  ();
void gfo ();
void gfo_adiab ();
void plm ();
void plm2 ();
void plm_adiab ();
void HydroKernel ();
void DustKernel ();
void DustDiffPresKernel ();
void DustDiffKernel ();
void Compute_Fluxes_Iso ();
void Compute_Fluxes_Adi ();
void Compute_Fluxes_pressureless1 ();
void Compute_Fluxes_pressureless2 ();
void Compute_Fluxes_Diffusion ();
real sgn();
real wavelimit();
void CorrectFluxesFromFinerLevel ();
void Source ();
void DustDiffusion ();
boolean boundary ();
real StoppingTimeLimit ();
real CourantLimit ();
real CourantLimitGlobal ();
void ViscousStress ();
void ApplyViscousStress ();
void prs_exit ();
void getgridsize ();
void getgridsizes ();
void WriteDescriptor ();
void BuildCommunicators ();
long StripCommunicators ();
void release_com ();
void DestroyCom ();
void DestroyHash ();
boolean CubeIntersect ();
void BlockButcher ();
boolean CubeIntersectPeriodic ();
void Comm_Alloc ();
void ExecCommSame ();
void ExecCommSameVar ();
void ExecCommUpVar ();
void ExecCommUpVarLIL ();
void ExecCommDownFlux ();
void ExecCommDownMean ();
void ExecCommUp ();
void ExecComm ();
void Write ();
void CheckBC ();
void TrueBC ();
void TrueBC_fp ();
void SetOverlapFlag ();
void ExecCommFlux ();
void MonitorConservative ();
void setfpe ();
void handfpe();
void GridBuild ();
void MakeDir ();
void ReadGrids ();
void CommAll ();
real Torque();
void TotalTorque();
void MonitorTorque ();
void MonitorTorqueZ ();
void keplerian_create();
real keplerian_init();
real keplerian_dust_init ();
real keplerian_dust_diffusion_init ();
real keplerian_gas_diffusion_init ();
real keplerian_pot();
real keplerian_comm();
real keplerian_boundary();
void GridAbs ();
void GridPos ();
void ReadField ();
void merge ();
FILE *prs_openrd();
FILE *prs_openr();
void merge_field ();
void setout ();
void ConstructGrids ();
void refine ();
void refine_field();
void FillSources ();
void FillSources_Dust ();
void FillSources_Predict ();
void FillSources_Predict_Dust ();
void FillSources_pot ();
void FillSources_geom_rad ();
void FillSources_geom_col ();
void FillSources2 ();
void FillSources_diff_dust ();
void FillSources_diff_gas ();
void memcpystride ();
void multarray ();
void arraymult ();
void FillSlopes ();
void FillSlopes_Dust ();
void mpiwarning ();
void AdjustBeamBoundaries ();
void addinitcode ();
void testdoubledefined ();
void testundef ();
void findcodeinit ();
void setinitcodes ();
void addpotcode ();
void testpotdoubledefined ();
void testundefpot ();
void findcodepot ();
void setpotcodes ();
void PressureCorrection ();
void SetCommSource ();
long nbvariables ();
long dimfield ();
void WriteCommSource ();
void WriteCommDest ();
PressureFaces *CreatePressureFaces ();
void FreePressureFaces ();
real InternalEnergy ();
real TotalInternalEnergy ();
void MonitorInternalEnergy ();
real KineticEnergy ();
real TotalKineticEnergy ();
void MonitorKineticEnergy ();
void MonitorTotalEnergy ();
real TorqueStockholm ();
void ApplyStockholmBoundaryConditions ();
void ApplyStockholmBoundaryConditionsDust ();
void DiskRadialDrift();
real **multiple_alloc_1D ();
void ReserveForBeam_PLM ();
void ReserveForSecondaryBeam_PLM ();
void ReserveForWorkPatch ();
void ReserveForSecondaryWorkPatch ();
void MakeBeam_PLM ();
void MakeSecondaryBeam_PLM ();
void AllocBeam_PLM ();
void AllocSecondaryBeam_PLM ();
void ReserveForBeam_MUSCL ();
void MakeBeam_MUSCL ();
void AllocBeam_MUSCL ();
void MakeCurrentFluidPatch ();
void MakeSecondaryFluidPatch ();
void Divergence ();
void FindBestSubcycling ();
void EvaluateBaseStepRatio ();
void EvaluateLevelCost ();
real SubCycleCost ();
void pInfo (const char *template, ...);
void pWarning (const char *template, ...);
void pError (const char *template, ...);
void ResetChronometer ();
void ReadChronometer ();
void FlushLog ();
void FlushInfo ();
void checkpara ();
int fselect ();
long GlobalInit ();
void ReadDefaultInOut ();
void SubsDef ();
void ExtractPath ();
FILE *FindGridFile ();
void DumpInitCode ();
void DumpPotCode ();
void ReadRedefined ();
void ParseRedefinedOptions ();
void muscl();
void Predictor ();
void convert_coord ();
void Ray_Integrate ();
void Ray_System_BoundingBox ();
void Ray_Tracing ();
void make_central_star ();
void make_planet ();
void make_background_stars ();
void AllocateSqueezedArrays ();
void CreateSqueezedArrays ();
void CreateAllSqueezedArrays ();
void HydroStatPrepare ();
SqueezedField *AllocateSqueezedArray ();
void WriteHydroStat ();
void WriteSqueezedField ();
void SetCenteredDensity ();
void SetInterfaceQuantities ();
void CorrectSourceTerm ();
void ForAllPatches ();
void SetToZero ();
void ResetPatch ();
real comm_adapt ();
real CorrKeplerianFlux ();
real IC_Average ();
real IC_2D_Mean ();
real IC_Corr ();
real IC_Density ();
real IC_Vrad ();
real IC_Vazimuth ();
real IC_Vcolatitude ();
real IC_Energy ();
real IC_Energy_tot ();
real GetInitialValue ();
real MassFluxRotating ();
real MassFluxSolidRotation ();
real MassFluxInertial ();
real EnergyTotFluxSolidRotation ();
real EnergyTotFluxInertial ();
real MomentumFluxRotating ();
real MomentumFluxInertial ();
real MomentumFluxSolidRotation ();
void cpTriplet ();
void resetTriplet ();
void swapl ();
void InitFluxCorrectionMeridian ();
real IC_Pressure ();
real IC_PressureByRad ();
boolean SquareIntersect ();
void GetMaxLevProjectedMeridian ();
void MultiFluid ();
void SetFluidProperties ();
void FluidCoupling ();
void MultifluidDiffusionPressure ();
void MultifluidDustEnergyToZero ();
void EnergyCorrection ();
void EnergyCorrection2 ();
void DustEnergyCorrection ();
real oplin();
real kappa();
void SubStep4 ();
void ComputeInitTemperatureField ();
void ComputeTemperatureField ();
void ComputeQplus ();
void AllocStellar();
void SubStep5 ();
void SolveMatrix();
void ComputeMatrixElements();
void ComputeEta();
void  ComputeStellarHeating();
void ComputeOpticalDepth();
void  ComputeOpacity();
void  ComputeOpacityInit();
void InitEnergyRad();
void reducegas();
void ExecCommOneVar ();
void ExecCommSameOneField ();
void GetOpticalDepthFromLevel (long);
void muscl_adiab ();
void Predictor_adiab ();
void Predictor_iso ();
void DensFloor ();
void DensFloorVelocityLimit ();
void DustDensFloor ();
Pair MinMax ();
void GetEnergyRadFromLevel ();
void RT_main ();
void PredictNewRadEnergy ();
void GetRadEnergyDifference ();
void MonitorSublimatedDustMass();
void MonitorRemovedDustMass();
void ReadRemDustMass();
void ReadSublDustMass();
boolean FileExists ();
real GetGamma();
void TestGammaSingleCell();
void TestGammaFluidPatch();
void WriteToFile();

