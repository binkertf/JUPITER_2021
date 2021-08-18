typedef double real;
typedef long   boolean;

struct param {
  char name[256];
  long type;
  char *variable;
  long read;
  long necessary;
};

typedef struct param Param;

struct hashparam {
  char name[MAXLINELENGTH];
  char stringvalue[MAXLINELENGTH];
  real floatvalue;
  long intvalue;
  long boolvalue;
};

typedef struct hashparam HashParam;

struct radixstruct {
	real	constant;
	real 	radius_slope;
	real	sintheta_slope;
};
typedef struct radixstruct RadixStruct;

enum {_density_, _energy_, _vrad_, _vazimuth_, _vcolatitude_, _other_};

enum {_Density_, _Energy_, _Velocity_, _vx2_, _vx3_,\
      _Potential_};
/* Note: if you add additional fields in the above enum, check straight away that
the Ptr[] array of FluidPatch is large enough (see gridtypes.h) */
