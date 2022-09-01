#define MAXLINELENGTH 512L
#define MAXNAMELENGTH 128L
#define YES 1L
#define TRUE 1L
#define NO 0L
#define FALSE 0L
#define REDEFINED 2L
#define REAL    1L
#define INT     0L
#define STRING  2L
#define BOOL    3L
#define MAXVARIABLES 1000L
#define MAXCHILDREN 100L
#define LOW 0L
#define HIGH 1L
#define INF 0L
#define SUP 1L
#define MIN 0L
#define MAX 1L
#define LEFT 0L
#define RIGHT 1L
#define CARTESIAN 0L		/* Do not change */
#define CYLINDRICAL 1L		/* these three */
#define SPHERICAL 2L		/* values. */
#define	PI	3.14159265358979323844
#define ONETHIRD 0.3333333333333333333
#define ARITHMETIC 0L
#define GEOMETRIC 1L
#define REFINEUPPERLIMIT 100L	/* Max refinement level ever. */
#define MASS 0L
#define MOMENTUM1LEFT 1L
#define MOMENTUM2LEFT 2L
#define MOMENTUM3LEFT 3L
#define ENERGY 4L
#define MOMENTUM1RIGHT 6L
#define MOMENTUM2RIGHT 7L
#define MOMENTUM3RIGHT 8L
#define MAXLENGTHONEDIM 512L
#define OUTSIDE 0L
#define INSIDE 1L
#define NGH 2L			/* Number of ghosts zones of meshes. 2 for Colella's method */
#define GHOST 0L
#define FLUX 1L
#define MEAN 2L
#define MAXGRIDS 1000l
#define PREDICT 0
#define UPDATE 1

#define SMALLESTDENSRATIO 1e-18 //For backward compatibility only

#define TEMPRATIO 1.0 // PLEASE TEST

/* Bitwise flag on communicator */
#define SRCGHOSTSAME 1
#define SRCGHOSTUP 2
#define SRCGHOSTSAMEEXTEND 4
#define SRCGHOSTUPEXTEND 8
#define DESTGHOSTSAME 1
#define DESTGHOSTUP 2
#define DESTGHOSTSAMEEXTEND 4
#define DESTGHOSTUPEXTEND 8
#define EVERYWHERE 16
#define NONE 0
#define AUTO 1
#define FULL 2
#define STANDARD 0L
#define EQUILIBRIUM 1L

// Definitions for radiative transfer:
#define       	XMSOL 		  1.989e+33  // Solar Mass in grams
#define       	XMH    		  1.67e-24  // Hydrogen Mass
#define       	AU      	  1.496e+13  // Astronomical Unit in cm 
#define       	BOLTZ   	  1.38e-16  // Boltzmann Constant
#define       	CLIGHT		  2.99792458e+10  // Light speed
#define       	GRAVC   	  6.67e-08  // Gravitation Constant
#define       	ARC        	  7.56e-15  // Radiation Constant
#define       	RGAS    	  8.314e+07  // Gas Constant
#define      	RSUN    	  (6.96e+10) /* Sun Radius in cm */  
#define         VOL0     	  (R0*R0*R0)
#define         XM0      	  (XMSOL * XMSTAR)
#define         RHO0     	  (XM0 / VOL0)
#define         TIME0    	  (sqrt(VOL0 / GRAVC / (XMSTAR*XMSOL)))
#define         V0      	  (R0 / TIME0)
#define         OPAC0   	  (R0*R0 / XM0)
#define         XNU0    	  (R0*V0)
#define         TEMP0    	  (V0*V0 / RGAS)
#define         P0       	  (RGAS * RHO0 * TEMP0)
#define         E0       	  (XM0 * R0*R0/ (TIME0*TIME0))
#define         SLUM     	  (1.0* TIME0 / E0)
#define    	C0	 	  (CLIGHT/V0)
#define     	A0	 	  (ARC*pow(TEMP0,4.0)/(RHO0*V0*V0))
#define      	SIGMARAD    	  (0.25*A0*C0)
#define        	BOUNDTEMP  	  3.0000000000000000000000000  // cooling temperature, we set it now equal to CMB temperature
#define        	MU        	  2.3000000000000000000000000 // mean molecular weight for standard solar mixture (Kley+ 2009)
#define        	CV      	  (1.000000000000000000000000000000/MU/(GAMMA-1.0)) // specific heat coeff (const. volume) code units (1.0 in the numerator instead of RGAS, latter would lead to CV in cgs units)
#define   	CPUOVERLAP   	  2	/* Zeus-like overlap kernel. 2:transport; 2: source, 1:viscous stress  */
#define    	TSTAR0   	 (TSTAR/TEMP0)
#define    	RSTAR0    	 (RSTAR*RSUN/R0)
#define    	FSTAR    	 (SIGMARAD*pow(TSTAR0,4.0)*RSTAR0*RSTAR0)
#define    	EPSILON    	 3.0e-6 // 5.e-6 how long the iteration should go in RT and Stellar
#define    	BoundEcode       (pow((BOUNDTEMP/TEMP0),4.0)*A0)


// Definitions for calculating GAMMA
#define CONST_h       6.62606876e-27     /**< Planck Constant.           */
#define CONST_eV      1.602176463158e-12 /**< Electron Volt in erg.      */
#define CONST_me      9.1093826e-28      /**< Electron mass.             */
#define CONST_amu     1.66053886e-24     /**< Atomic mass unit.          */
#define CONST_mp 1.67262171e-24 /**< Proton mass. */
#define H_MASS_FRAC       0.7110
#define  He_MASS_FRAC  (1 - H_MASS_FRAC) /* Effective Y and not 0.2741*/
/* Vibrational and Rotational Temperatures for molecular hydrogen  */
#define THETA_V  6140.0
#define THETA_R    85.5
#define T_CUT_RHOE  10.0
#define DEG_x      0  /* hydrogen ionization degree */
#define DEG_y      1  /* molecular hydrogen dissociation degree */
#define DEG_z1     2  /* helium single ionization degree */
#define DEG_z2     3  /* helium double ionization degree */
#define HELIUM_IONIZATION   NO
#define H_DISSOCIATION NO
#define H_IONIZATION NO
#define TESTING_GAMMA NO

/*! Define the conversion constant between dimensionless
    temperature prs/rho and physical temperature in Kelvin,
    T = (prs/rho)*KELVIN*mu
*/
#define KELVIN (V0*V0*CONST_amu/BOLTZ)


//#define ROTFRAME sqrt(1.0-2./7.*0.05*0.05)
#define ROTFRAME 0.0


#define INSPECT_REAL( var) {printf ("%s = %.18g\n", #var, var);}
#define INSPECT_INT( var) {printf ("%s = %d\n", #var, var);}

#define CHECK( m) {real temp; \
temp = CurrentFluidPatch->Velocity[1][m];\
printf ("Rad velocity in zone %d: %g at file %s and line %d\n", m, temp, __FILE__, __LINE__);\
if (fabs(temp) > 1e-12) exit(1);\
}

#define JUP_SAFE( statement)  statement;

#define DUMP( string, nb) WriteWorkArray (CurrentFluidPatch->string, #string, nb);

/*
long mmm;
long allj;
long currentoutnb;

#define JUP_SAFE( statement)  {		\
statement;	\
printf ("Now at %d, executing %s\n", currentoutnb, #statement);			\
for (mmm = 0 ; mmm < (CurrentFluidPatch->desc->gncell[0])*(CurrentFluidPatch->desc->gncell[1])*(CurrentFluidPatch->desc->gncell[2]);mmm++) {\
if (isnan(CurrentFluidPatch->Density[mmm])) {\
printf ("Density is Nan at m=%ld after call to %s in file %s at line %d\n", mmm, #statement, __FILE__, __LINE__); \
allj=1;\
}\
if (isnan(CurrentFluidPatch->Energy[mmm])) {\
printf ("Energy is Nan at m=%ld after call to %s in file %s at line %d\n", mmm, #statement, __FILE__, __LINE__); \
allj=1;\
}\
if (isnan(CurrentFluidPatch->Velocity[0][mmm])) {\
printf ("Vel0 is Nan at m=%ld after call to %s in file %s at line %d\n", mmm, #statement, __FILE__, __LINE__); \
allj=1;\
}\
if (isnan(CurrentFluidPatch->Velocity[1][mmm])) {\
printf ("Vel1 is Nan at m=%ld after call to %s in file %s at line %d\n", mmm, #statement, __FILE__, __LINE__); \
allj=1;\
}\
if (isnan(CurrentFluidPatch->Velocity[2][mmm])) {\
printf ("Vel2 is Nan at m=%ld after call to %s in file %s at line %d\n", mmm, #statement, __FILE__, __LINE__); \
allj=1;\
}\
if (allj ==1) exit (1);\
}\
}*/

Pair resmm;
#define MINMAX( array) {\
resmm = MinMax (array);\
 pInfo ("Array %s: min/max = %g / %g\n", #array, resmm.min, resmm.max);	\
}\

#define CHECK_CONVERGENCE

//#define CHECK_CONVERGENCE {\
//if ((prs < 0.0) || (isnan(prs))) {\
//  fprintf (stderr, "CONVERGENCE PROBLEM IN %s AT LINE %d\n", __FILE__, __LINE__);\
//  INSPECT_REAL (rhoL);\
//  INSPECT_REAL (rhoR);\
//  INSPECT_REAL (uL);\
//  INSPECT_REAL (uR);\
//  INSPECT_REAL (aL);\
//  INSPECT_REAL (aR);\
//  INSPECT_REAL (prs);\
//  exit (EXIT_FAILURE);\
// }\
//}
