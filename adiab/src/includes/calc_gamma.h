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


#define HELIUM_IONIZATION   NO

#define DEG_x      0  /* hydrogen ionization degree */
#define DEG_y      1  /* molecular hydrogen dissociation degree */
#define DEG_z1     2  /* helium single ionization degree */
#define DEG_z2     3  /* helium double ionization degree */

#define TEMP      0  /* temperature */
#define RHO      1  /* density */

double Gamma1(double *v);
