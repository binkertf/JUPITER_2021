#include "jupiter.h"

/* In this function we adapt the value that is sent via communicators
(to other levels, generally, be it via a Ghost or Mean kind of
communication). In the case in which hydrostatic equilibrium is
enforced, the perturbed density rather than the density is
communicated. If furthermore KEPLERIAN is set, we communicate
Omega*r^1.5 rather than Omega.  In the case in which Hydrostatic
Equilibrium is not enforced, we communicate quantities according to
the transformations defined in keplerian.c::keplerian_comm. A
direction of +1 is required to fill the communicator before a
communication is initiated, while a direction of -1 restores the
normal field values from a communicator upon communication
completion. Passing both i[3] and m in the arguments is redundant, but
both are usually available from the calling function, so they are both
sent for optimization purposes. */

inline real comm_adapt (value, component, m, i, sqz, desc, dir, fluidrank)
     real value;		/* value to be adapted */
     long component;		/* type of this value (density, etc.) */
     long m;			/* Location on grid (1D index) */
     long i[3];			/* Coordinates on grid */
     long sqz[3]; 		/* Squeezed stride */
     tGrid_CPU *desc;		/* Associated grid descriptor */
     int dir;	       /* direction +1: before comm, -1: after comm */
     long fluidrank;
{
  long ms;
  real rhoc, enec, result;
  real *rad, *col;
  FluidPatch *fluid;
  result = value;
  /* We have to make the distinction between 'HydroStaticEnforce',
     which is set as soon as the parameters are read, and
     'HydroStaticReady', which is set once the squeezed arrays
     describing the equilibrium are available. The reason for that is
     that initially all the meshes are synchronized by communications,
     at the very begining of the hydrostatic initialization process,
     when the squeezed arrays are not created yet.
 */
  fluid = desc->Fluid;
  while (fluidrank > fluid->FluidRank) fluid=fluid->next;
  if (HydroStaticReady) {
    if (component == _density_) {
      ms = i[0]*sqz[0]+i[1]*sqz[1]+i[2]*sqz[2];
      rhoc = fluid->Rho_eq_c->field[ms];
      result = value - (real)dir*rhoc;
    }
    if ((component == _energy_) && (!Isothermal)) {
      ms = i[0]*sqz[0]+i[1]*sqz[1]+i[2]*sqz[2];
      enec = fluid->Ene_eq_c->field[ms];
      result = value - (real)dir*enec;
    }
  }
  if ((KEPLERIAN) && ((!HydroStaticReady) || (component == _vazimuth_))) {
    rad = desc->Center[_RAD_];
    col = desc->Center[_COLAT_];
    result = keplerian_comm(value, component, rad[m], col[m], dir);
  }
  return result;
}
