/** The mesh structure that contains the information independent on
    the number of processing elements */

struct tgrid {
  long level;			/**< Refinement level */
  long ncell[3];		/**< Number of zones in each dimension (excluding ghosts) */
  long gncell[3];		/**< Number of zones in each dimension (including ghosts) */
  long nsize[3];		/**< Absolute size (excluding ghosts) */
  long ncorner_min[3];		/**< Absolute position (min corner, excluding ghosts) */
  long ncorner_max[3];		/**< Absolute position (max corner, excluding ghosts) */
  long gnsize[3];		/**< Absolute size (including ghosts) */
  long gncorner_min[3];		/**< Absolute position (min corner, including ghosts) */
  long gncorner_max[3];		/**< Absolute position (max corner, including ghosts) */
  real size[3];			/**< Size in physical units */
  real corner_min[3];		/**< Position of min corner in physical units */
  real corner_max[3];		/**< Position of max corner in physical units */
  real *Edges[3];		/**< Interface position in each dimension */
  boolean present;		/**< True if and only if that grid is mapped */
				/**< by local CPU */
  long number;			/**< Logical number of grid */
  struct tgrid *next;		/**< Pointer to next grid */
  struct tgrid *prev;		/**< Pointer to previous grid */
  long BoundaryConditions[3][2];/**< List of BC's on faces */
  long minCPU;			/**< range of processes over which */
  long maxCPU;			/**< the grid is distributed */
  long dn[3];			/**< Cell size in absolute units */
  long linenumber;		/**< Line number in grid file */
  long monoCPU;			/**< Information used for a restart only */
};

typedef struct tgrid tGrid;
typedef struct fluidpatch FluidPatch;

/** The information relative to a fluid patch that is accessible to a
given processing element.  The 'parent' (or 'Parent') information
refers to the tGrid of which the current FluidPatch is the local
subset. */

struct tgrid_cpu {
  long   parent;       		/**< Logical number of "parent" grid */
  long   number;		/**< Internal logical number */
  tGrid *Parent;		/**< Pointer to "parent */
  long   level;			/**< Refinement level */
  int    cpu;			/**< process number of corresponding CPU grid */
  long   ncell[3];		/**< Number of zones in each dimension (excluding ghosts)*/
  long   gncell[3];		/**< Number of zones in each dimension (including ghosts)*/
  long   stride[3];		/**< Data storage stride in each dimension */
  long   nsize[3];		/**< Absolute size */
  long iface[3][2];		/**< Interface type (CPU->brother, True BC, Mesh->other mesh) */
  long ncorner_min[3];		/**< Absolute position (min corner) */
  long ncorner_max[3];		/**< Absolute position (max corner) */
  long pcorner_min[3];		/**< Position in parent (in cells) - min */
  long pcorner_max[3];		/**< Position in parent (in cells) - max */
  long gnsize[3];		/**< Absolute size (including ghosts) */
  long gncorner_min[3];		/**< Absolute position (min corner, including ghosts) */
  long gncorner_max[3];		/**< Absolute position (max corner, including ghosts) */
  real size[3];			/**< Size in physical units */
  real corner_min[3];		/**< Position of min corner in physical units */
  real corner_max[3];		/**< Position of max corner in physical units */
  real *InvVolume;		/**< Inverse of volume of each cell */
  real *Center[3];		/**< Center of zones in each dimension */
  real *Metric[3][2];		/**< 1D Arrays for metric coefficients */
  real *InvMetric[3][2];	/**< 1D Arrays for metric coefficients */
  real *InterSurface[3];	/**< Surface of interfaces in each dimension */
  real *Edges[3];		/**< Interface position in each dimension */
  char *Hidden;			/**< A set of flags which says if the zone lies behind another */
  char *CommSource;		/**< A set of bit flags used by LeVeque method */
  char *CommDest;		/**< A set of bit flags used by LeVeque method */
  char *SrcCommPerp;		/**< A set of bit flags to switch between LeVeque/Keplerian */
  char *DestCommPerp;		/**< A set of bit flags to switch between LeVeque/Keplerian */
  long *MaxLevMeridianProj;	/**< Array of azimuthally projected highest level (see merid*\.c for details)  */
  FluidPatch *Fluid;		/**< Pointer to the first fluid patch (follow with ->next for multifluid) */
  int color;                    /**< Variables required for the Stellar Irradiation with MPI */
  int colorz;                   /**< Variables required for the Stellar Irradiation with MPI */
  int key;                      /**< Variables required for the Stellar Irradiation with MPI */
  real *optical_depth;          /**< Total optical depth of the patch at stellar wavelengths along radial direction */
  struct tgrid_cpu *next;	/**< Pointer to next CPU Grid */
  struct tgrid_cpu *prev;	/**< Pointer to previous CPU Grid */
};

typedef struct tgrid_cpu tGrid_CPU;

/** Container used to synchronize the ghost zones between different
grids.  A given communicator is used by two different CPU_Grids (and
two only), and it contains all the information needed to pick up
values in the source CPU_grid an affect them in the destination
CPU_grid. It handles all the hydrodynamics variables at once, in order
to lower the communication cost.
 */

struct communicator {
  int tag;			/**< Tag used by MPI. Necessarily
				   'int' and not 'long' as this has to
				   correspond to MPI implementation */
  long imin_dest[3]; /**< lower boundary of the communicator on the
				     destination CPU_Grid. The lower
				     boundary is included.
		     */
  long imax_dest[3]; /**< upper boundary of the communicator on the
				     destination CPU_Grid. The upper
				     boundary is excluded.
		     */
  long imin_src[3]; /**< lower boundary of the communicator on the
				     source CPU_Grid. The lower
				     boundary is included.
		    */
  long imax_src[3]; /**< upper boundary of the communicator on the
				     source CPU_Grid. The lower
				     boundary is included.
		    */
  long dest_level; /**< Refinement level of the destination grid */
  long src_level; /**< Refinement level of the source grid */
  long grid_src;  /**< Logical number of the parent of the source grid */
  long grid_dest; /**< Logical number of the parent of the destination grid */
  long facedim; /**< Dimension perpendicular to the communicator */
  long faceside; /**< INF or SUP, relative to the destination */
  long CPU_src, CPU_dest; /**< processing element number of source and destination */
  long nb_src, nb_dest; /**< Logical numbers of the source and destination grids */
  long Imin[3], Imax[3]; /**< Corners of the communicator in absolute coordinates */
  long type;			/**< Should be one of these values: GHOST, MEAN, FLUX */
  real *buffer; /**< Container for the data to be transferred */
  tGrid_CPU *srcg, *destg; /**< Pointers to the descriptors of the
			      source and destination CPU_Grids */
  long size; /**< size of communicator, for the destination mesh. The
		size of a communicator is always its size for the
		destination mesh. It is expressed in number of
		zones. */
  struct communicator *next; /**< Pointer to the next communicator */
  struct communicator *prev; /**< Pointer to the previous communicator */
};

/** Hash tables of communicators are built to find efficiently
redundant communicators. A given hash table groups all the
communicators that have the same source CPU_Grid. The tables
are handled as chained list. */

struct commhash {
  struct communicator *com; /**< Communicator of the current hash table element */
  struct commhash *next;    /**< Pointer to the next element of the hash table */
  struct commhash *prev;    /**< Pointer to the previous element of the hash table */
};

typedef struct communicator Communicator;

typedef struct commhash CommHash;

/** Structure that handles scalar hydrodynamics fields (eg density) */

struct scalarfield {
  char *Name;      /**< Name of the field (eg "density", "energy") */
  real *Field;     /**< Pointer to the array that contains the field */
  tGrid_CPU *desc; /**< Pointer to the associated descriptor of the CPU_Grid */
};

typedef struct scalarfield ScalarField;

/** Structure that handles vectorial hydrodynamics fields (eg velocity) */

struct vectorfield {
  char *Name;       /**< Name of the field (eg "velocity") */
  real *Field[3];   /**< Pointer to the array that contains the field */
  tGrid_CPU *desc;  /**< Pointer to the associated descriptor of the CPU_Grid */
};

typedef struct vectorfield VectorField;

/** Structure that is used to monitor fluxes across the faces of a given CPU_Grid */

struct interfaceflux {
  real *Flux[3][2]; /**< Pointer to the 6 arrays of fluxes (3 dimensions * 2 sides) */
  tGrid_CPU *desc;  /**< Pointer to the associated descriptor */
};

typedef struct interfaceflux InterfaceFlux;

/** Structure that is used to monitor the pressure on the faces of a given CPU_Grid */

struct pressurefaces {
  real *Pressure[3][2];  /**< Pointer to the 6 arrays of pressure (3 dimensions * 2 sides) */
  tGrid_CPU *desc;	 /**< Pointer to the associated descriptor */
};

typedef struct pressurefaces PressureFaces;

/** Structure that is used to store specific configurations which can have symmetry properties.
Such configurations are essentially equilibrium configurations. In that case the symmetry of the
system can be used to reduce the amount of data that needs to be stored: the array are "squeezed"
along the invariant dimension (e.g. the azimuth is an invariant dimension for a Keplerian disk
in hydrostatic equilibrium) */

struct squeezedfield {
  char *Name;			/**< Name of the field  */
  tGrid_CPU *desc;		/**< Pointer to the associated descriptor */
  long gncell[3];		/**< Size of the squeezed field. It is 1 along a squeezed dimension. */
  long stride[3];		/**< Associated stride */
  real *field;			/**< Pointer to the array that contains the field */
};

typedef struct squeezedfield SqueezedField;

/** Structure that is used to store the information relative to the hydrodynamics
    field of all CPU_Grids. */

struct fluidpatch {
  char *Name;			/**< The name of the fluid (eg gas, or dust, etc) */
  ScalarField *Density;         /**< Pointer to the field \f$\rho\f$*/
  real *StartField;             /**< Pointer to the beginning of the first array (they are all contiguous in memory) */
  VectorField *Velocity;        /**< Pointer to the velocity field */
  ScalarField *Energy;          /**< Pointer to the volumic energy field */
  ScalarField *TotalEnergy;     /**< Pointer to the volumic total energy field */
  ScalarField *OpticalDepth;    /**< Pointer to the optical depth integrated along the star direction */
  ScalarField *EnergyRad;       /**< Pointer to the radiative energy field */
  ScalarField *StellarHeating;  /**< Pointer to something obvious */
  ScalarField *EradDeriv;       /**< Time derivative of radiative energy */
  ScalarField *Temperature;       /**< Temperature */
  ScalarField *Gamma;           /**< Gamma */
  ScalarField *Opacity;       /**< Opacity */
  ScalarField *TauCell;       /**< Optical Depth value in a given cell */
  ScalarField *Diffcoeff;        /**< DiffCoefficients */
  real *Ptr[23];		/**< Short cut to fields */
  InterfaceFlux  *MassFlux;     /**< Pointer to the structure that handles the mass flux at the grid boundary */
  InterfaceFlux  *DiffFlux;
  InterfaceFlux  *MomentumFlux[3]; /**< Pointer to the 3 structures that handles the momenta fluxes at the grid boundary */
  InterfaceFlux  *TotalEnergyFlux;   /**< Pointer to the structure that handles the total energy flux at the grid boundary */
  InterfaceFlux  *EnergyFlux;   /**< Pointer to the structure that handles the internal energy flux at the grid boundary */
  PressureFaces  *Pressure;     /**< Pointer to the structure that handles the arrays of pressure at the grid boundary */
  ScalarField *Potential;	/**< Pointer to the scalar field
				   potential. It is stored so as to
				   avoid an evaluation at each time
				   step, in the case it is constant in
				   time.  */
  SqueezedField *Rho_eq_c;	/**< Field that handles the zone centered density at hydrostatic equilibrium */
  SqueezedField *Ene_eq_c;	/**< Field that handles the zone centered internal energy at hydrostatic equilibrium */
  SqueezedField *Rho_eq_i[3];	/**< 3 fields that handle the interface density at hydrostatic equilibrium */
  SqueezedField *Ene_eq_i[3];	/**< 3 fields that handle the interface internal energy at hydrostatic equilibrium */
  SqueezedField *Cs2_i[3]; 	/**< 3 fields that handle the interface energy at hydrostatic equilibrium */
  SqueezedField *Pres_eq_i[3];	/**< 3 fields that handle the interface pressure at hydrostatic equilibrium */
  SqueezedField *Source_corr[3];/**< 3 fields that contain the additive source term correction at hydrostatic equilibrium */
  tGrid_CPU      *desc;		/**< Pointer to the descriptor associated to this fluid patch  */
  boolean PotentialSet;		/**< Indicates whether the potential has already been evaluated  */
  real *MassFluxCorrection1;	/**< Array 1 of azimuthal mass flux correction, if needed */
  real *MassFluxCorrection2;    /**< Array 2 of azimuthal mass flux correction, if needed */
  real *MomentumFluxCorrection1;	/**<  Idem for the angular momentum */
  real *MomentumFluxCorrection2;	/**< See merid*\.c for details, as well as sample.c */
  real *EnergyFluxCorrection1;		/**< (look for C_mass_1, C_mass_2, C_mom_1 and C_mom_2) */
  real *EnergyFluxCorrection2;
  FluidPatch *next;		/**< For multifluid calculations, pointer to the next fluid  */
  FluidPatch *prev;		/**< For multifluid calculations, pointer to the previous fluid  */
  long FluidRank; 		/**< Fluid number for multifluid calculation  */
  long InitCode;		/**< Initialization code for this fluid  */
  long InitCodeEq;		/**< Initialization code for this fluid (equilibrium configuration) */
  boolean PreviousEradExists;
};

struct fluidwork {
  real *Density;
  real *Density_Pred;
  real *Energy_Pred;
  real *Potential;
  real *Energy;
  real *Energy_tot;
  real *Divergence;
  real *Momentum[3];
  real *Velocity[3];
  real *Velocity_Pred[3];
  real *SourceDiv;		/* to store geometric term */
  real *SourceRhoPred;		/* to store geometric term */
  real *Accel[3];
  real *SourceVelocity[3];
  real *SlopeDensity[3];
  real *SlopeEnergy[3];
  real *SlopePotential[3];
  real *SlopeVelocity[3][3];
  real *InterfaceVel[3][3]; 	/* So memory expensive... */
  real *StressTensor[3];
  real *Flux_mass[3];
  real *Flux_diff[3];
  real *Flux_energy[3];
  real *Centrifugal_acc;
  real *Coriolis;
  real *RawMassFlux[3];
  real *Flux_tot_energy[3];
  real *Flux_mom[3][3];		/* right index: dim for face, left index: comp of mom */
  real *InterfacePressure[3];		/* Pressure at interfaces */
  real *Temperature;                /* From now on list arrays which are used in Radiative Transfer */
  real *Gamma;
  real *Qplus;                              /* Qplus: viscous heating */
  real *Diffcoeff;                          /* Aar...Rhs are coefficents of the matrix of the Eq. B.6*/
  real *Aarmatrix;                          /*in Kley et al. 2009, Appendix */
  real *Aatmatrix;
  real *Aapmatrix;
  real *Ccrmatrix;
  real *Cctmatrix;
  real *Ccpmatrix;
  real *Bbmatrix;
  real *Rhsmatrix;
  real *TemperOld;
  real *TemperNew;
  real *EnergyNew;
  real *DiffTemp;
  real *OpaR;
  real *OpaP;
  real *OpaS;
  real *StellarRad;
  real *Eta1;
  real *Eta2;
  real *EnergyRad;
  real *EnergyRadNew;
  real *EnradOld;
  real *DiffEnrad;
  real *TauOptical;
  real *TauMax;
  real *GlobalTauMax;
  real *TauCell;
  tGrid_CPU *desc;
  FluidPatch *Fluid;
  FluidPatch *next;		/* For multifluid calculations, pointer to the next fluid  */
  long Size[3];
  long GroupSize[3];
};

typedef struct fluidwork FluidWork;

struct beam {
  real *rho;
  real *u;
  real *v_perp[2];
  real *cs;
  real *slope_rho[3];
  real *slope_energy[3];
  real *slope_u[3];
  real *slope_v_perp[2][3];
  real *center;
  real *rawcoord;
  real *edge;
  real *srcrho;
  real *intersurface;
  real *mass_flux;
  real *diff_flux;
  real *energy_flux;
  real *tot_energy_flux;
  real *momentum_flux[3];
  real *rhoL;
  real *rhoR;
  real *uL;
  real *uR;
  real *u_interface[3];
  real *eL;
  real *eR;
  real *rho_pred;
  real *e_pred;
  real *u_pred;
  real *v_perp_pred[2];
  real *v_perp_L[2];
  real *v_perp_R[2];
  real *pressure_godunov;
  real *metperp[2];
  real *source;
  real *source_perp[2];
  real *HS_cent_rho;
  real *HS_int_rho;
  real *HS_cent_ene;
  real *HS_int_ene;
  real *cs2i;
  real *rad_arr;
  real *sin_theta_arr;
  real *centrifugal_acc;
  real *coriolis;
  real *raw_mass_flux;
  long dim[3];
  long length;
  tGrid_CPU *desc;
  boolean true_bc[2];
  boolean ZoneSplit;
  real rawcoord1;
  real rawcoord2;
  real radius;
  real colatitude;
  real radmin;
  real radmax;
  real masscorr1;
  real masscorr2;
  real momcorr1;
  real momcorr2;
  real enercorr1;
  real enercorr2;
};

typedef struct beam Beam;

struct gridfileinfo {
  long number;
  long linenumber;
  long nc_min[3];
  long nc_max[3];
  long size[3];
  real xmin[3], xmax[3];
  long level;
  long bc[6];
  long monoCPU;			/* Restart information */
  boolean last;
};

typedef struct gridfileinfo GridFileInfo;
