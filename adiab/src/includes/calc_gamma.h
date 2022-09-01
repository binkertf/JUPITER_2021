#define CONST_h       6.62606876e-27     /**< Planck Constant.           */
#define CONST_eV      1.602176463158e-12 /**< Electron Volt in erg.      */
#define CONST_me      9.1093826e-28      /**< Electron mass.             */
#define CONST_amu     1.66053886e-24     /**< Atomic mass unit.          */

#define H_MASS_FRAC       0.7110
#define  He_MASS_FRAC  (1 - H_MASS_FRAC) /* Effective Y and not 0.2741*/

/*! Define the conversion constant between dimensionless
    temperature prs/rho and physical temperature in Kelvin,
    T = (prs/rho)*KELVIN*mu
*/
#define KELVIN (V0*V0*CONST_amu/BOLTZ)

#define DEG_x      0  /* hydrogen ionization degree */
#define DEG_y      1  /* molecular hydrogen dissociation degree */
#define DEG_z1     2  /* helium single ionization degree */
#define DEG_z2     3  /* helium double ionization degree */

#define HELIUM_IONIZATION   NO
#define H_DISSOCIATION NO
#define H_IONIZATION NO

/* Vibrational and Rotational Temperatures for molecular hydrogen in Kelvin  */
#define THETA_V  6140.0 /* 0.529104 eV */
#define THETA_R    85.5 /* 0.0073678 eV */

#define T_CUT_RHOE  10.0

#define NONZERO_INITIALIZE YES /* Fill arrays to nonsense values to catch uninitialized values later in the code */
char     *Array1D (int, size_t);
#define ARRAY_1D(nx,type)          (type    *)Array1D(nx,sizeof(type))

/* Switch to choose the molecular hydrogen spin states.
 * 0 : Only Para hydrogen ,
 * 1 : Equilibrium (DEFAULT),
 * 2 : Ortho to para ratio of 3:1.                                 */
#define ORTHO_PARA_MODE 2


double Gamma1(double temperature, double pressure);

