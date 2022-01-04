/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varsubst" for details       */
/*                              */
/********************************/
extern boolean Verbose, Merge, AddSubPatch;
extern boolean Disable, Restart, CFLMonitor, Refine[3];
extern boolean WriteGhosts, Isothermal, MonitorCons;
extern boolean SuperImpose, AllCPUs, MonitorTimeSpent;
extern boolean TorqueMon, EneMon, NoStockholm, TorqStock;
extern boolean AllowFlushLog, QuietOutput, RedefineOptions;
extern boolean RayTracing, RayMirror, CoarseRayTracing;
extern boolean HydroStaticEnforce, HydroStaticReady;
extern boolean Radiative, HalfDisk, Stellar;
extern boolean Stretch, SmoothTaper, constSt;
extern long DebugOutputNb;
extern boolean mPLM, mGFO, mMUSCL;
extern boolean constDTG;
extern boolean __SPHERICAL, __CARTESIAN, __CYLINDRICAL;
extern boolean FreshStart;
extern long SpacingDim[3], InvCoordNb[3], CoordNb[3];
extern long SubCycling;
extern long NbFluids;
extern char FluidName[10][MAXLINELENGTH], InitCodeNames[10][MAXLINELENGTH], InitCodeNamesEq[10][MAXLINELENGTH];
extern long CoordType;
extern long InitMode;
extern long _Radial_, _Azimuthal_, _Vertical_, _Colatitude_;
extern long _RAD_, _AZIM_, _COLAT_;
extern long _X_, _Y_, _Z_, LevMax;
extern void    (*__Riemann_Solver)();
extern void    (*__Prepare_Riemann_States)();
extern void    (*__Compute_Fluxes)();
extern void    (*__Compute_Fluxes_pressureless)();
extern void    (*__ExecCommUpVar)();
extern real     GlobalDate;
extern int CPU_Rank, CPU_Number;
extern tGrid_CPU *Grid_CPU_list;
extern long Ncell0[3];
extern real corner_min0[3], corner_max0[3];
extern long ncorner_min0[3], ncorner_max0[3];
extern long Nghost[3], InitCode, PotentialCode;
extern boolean Periodic[3];
extern tGrid *GridList;
extern Communicator *ComListGhost;
extern Communicator *ComListFlux;
extern Communicator *ComListMean;
extern CommHash **CommHashSrc;	     /* This hash table is used for ghost com. only */
extern long CurrentOutputNumber;
extern FluidWork *CurrentFluidPatch;
extern FluidWork *SecondaryFluidPatch;
extern real _Small_Rho;
extern real _Small_Energy;
extern real TOGPO, TOGMO, GMOGPO, OOG, OOGMO;
extern char OutputDir[MAXLINELENGTH];
extern char DefaultOut[MAXLINELENGTH];
extern char DefaultIn[MAXLINELENGTH];
extern char InputPath[MAXLINELENGTH];
extern char ShortParName[MAXLINELENGTH];
extern char SubPatchInfo[MAXLINELENGTH];
extern char RayTracingInfo[MAXLINELENGTH];
extern boolean MomentumCorrection[3][3];
extern real GridFriction[3];
extern real DatePotentialConstant;
extern real GlobalDateInit;
extern boolean EmbeddedGridFile;
extern char EmbeddedGridFileName[MAXLINELENGTH];
extern long MaxLevel;
extern long NDimHydroStat;		/* Number of dimensions along which
extern 				   hydrostatic equilibrium is
extern 				   enforced */
extern long CorrHydroStat[3];		/* Is 1 for each dimension along which
extern 				   hydrostatic equilibrium is
extern 				   enforced, and 0 otherwise */
extern long dimcorr[3];		/* Gives the coordinate number of the
extern 				   NDimHydroStat coordinates along
extern 				   which hydrostatic equilibrium is
extern 				   enforced */
extern long SqueezeDim[3];		/* Is 1 if system in equilibrium is
extern 				   invariant along dimension (0,1 or
extern 				   to reduce memory waste by the
extern 				   arrays of hydrostatic corrections,
extern 				   as well as cache exceptions that
extern 				   would ensue */
extern long InitCode_Eq;		/* Initial conditions code of
extern 				   hydrostatic equilibrium */
extern long PotCode_Eq;		/* Potential code of
extern 				   hydrostatic equilibrium */
