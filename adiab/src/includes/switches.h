boolean Verbose = NO, Merge = NO, AddSubPatch = NO;
boolean Disable = NO, Restart = NO, CFLMonitor = NO, Refine[3];
boolean WriteGhosts = NO, Isothermal = NO, MonitorCons = YES;
boolean SuperImpose = NO, AllCPUs = NO, MonitorTimeSpent = NO;
boolean TorqueMon = NO, EneMon = NO, NoStockholm = NO, TorqStock = NO;
boolean AllowFlushLog = NO, QuietOutput = NO, RedefineOptions = NO;
boolean RayTracing = NO, RayMirror = NO, CoarseRayTracing = NO;
boolean HydroStaticEnforce = NO, HydroStaticReady = NO;
boolean Radiative = NO, HalfDisk = NO, Stellar = NO;
boolean Stretch = NO, SmoothTaper = NO, constSt = NO;
long DebugOutputNb = 0;
/* See comm_adapt.c for the use of HydroStaticReady */
boolean mPLM = NO, mGFO = NO, mMUSCL = NO;
boolean constDTG = NO;
boolean __SPHERICAL = NO, __CARTESIAN = NO, __CYLINDRICAL = NO;
boolean FreshStart = YES;
long SpacingDim[3], InvCoordNb[3] = {0, 1, 2}, CoordNb[3] = {0, 1, 2};
long SubCycling;
long NbFluids;
char FluidName[10][MAXLINELENGTH], InitCodeNames[10][MAXLINELENGTH], InitCodeNamesEq[10][MAXLINELENGTH];
long CoordType = CARTESIAN;
long InitMode = STANDARD;
long _Radial_, _Azimuthal_, _Vertical_, _Colatitude_;
long _RAD_, _AZIM_, _COLAT_;
long _X_, _Y_, _Z_, LevMax=0;
void    (*__Riemann_Solver)();
void    (*__Prepare_Riemann_States)();
void    (*__Prepare_Riemann_States_Dust)();
void    (*__Compute_Fluxes)();
void    (*__Compute_Fluxes_pressureless)();
void    (*__Compute_Fluxes_Dust)();
void    (*__ExecCommUpVar)();
real    (*__TVDslope)();
real     GlobalDate = 0.0;
int CPU_Rank, CPU_Number;
tGrid_CPU *Grid_CPU_list = NULL;
/* Template grid definition */
long Ncell0[3];
real corner_min0[3], corner_max0[3];
long ncorner_min0[3], ncorner_max0[3];
long Nghost[3], InitCode, PotentialCode;
boolean Periodic[3];
tGrid *GridList = NULL;
Communicator *ComListGhost = NULL;
Communicator *ComListFlux = NULL;
Communicator *ComListMean = NULL;
CommHash **CommHashSrc;	     /* This hash table is used for ghost com. only */
long CurrentOutputNumber = 0;
FluidWork *CurrentFluidPatch = NULL;
FluidWork *SecondaryFluidPatch = NULL;
real _Small_Rho = 1e20;
real _Small_Energy;

real TOGPO, TOGMO, GMOGPO, OOG, OOGMO;

char OutputDir[MAXLINELENGTH];
char DefaultOut[MAXLINELENGTH];
char DefaultIn[MAXLINELENGTH];
char InputPath[MAXLINELENGTH];
char ShortParName[MAXLINELENGTH];
char SubPatchInfo[MAXLINELENGTH];
char RayTracingInfo[MAXLINELENGTH];
boolean MomentumCorrection[3][3];
real GridFriction[3];
real DatePotentialConstant;
real GlobalDateInit;
boolean EmbeddedGridFile = NO;
char EmbeddedGridFileName[MAXLINELENGTH];
long MaxLevel=100L;
long NDimHydroStat;		/* Number of dimensions along which
				   hydrostatic equilibrium is
				   enforced */
long CorrHydroStat[3];		/* Is 1 for each dimension along which
				   hydrostatic equilibrium is
				   enforced, and 0 otherwise */
long dimcorr[3];		/* Gives the coordinate number of the
				   NDimHydroStat coordinates along
				   which hydrostatic equilibrium is
				   enforced */
long SqueezeDim[3];		/* Is 1 if system in equilibrium is
				   invariant along dimension (0,1 or
				   2), and 0 otherwise. This is used
				   to reduce memory waste by the
				   arrays of hydrostatic corrections,
				   as well as cache exceptions that
				   would ensue */
long InitCode_Eq;		/* Initial conditions code of
				   hydrostatic equilibrium */
long PotCode_Eq;		/* Potential code of
				   hydrostatic equilibrium */
