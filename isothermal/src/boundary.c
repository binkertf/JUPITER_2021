#include "jupiter.h"

inline boolean boundary (rho,e,u,v,w,rhog,eg,ug,vg,wg,x,xg,yg,zg,condition,line,predictive)
     real rho,e,u,v,w,x,xg,yg,zg;
     boolean predictive;	/* Predictive is true whenever
     'boundary' is called from within a Riemann problem predictor. It
     is false otherwise (ie when filling the mesh boundaries) */
     real *rhog,*eg,*ug,*vg,*wg;
				/* x and xg are the zone coordinate
                                   along dimension perpendicular to
                                   the boundary, respectively for
                                   source and dest (ie ghost) zones */
     long condition,line;		/* u is the perpendicular velocity */
				/* v and w are the transverse velocities */
{
  boolean Sym = TRUE;
  switch (condition) {
  /* usual conditions */
  case 1:	       /* REFLECTING */
    *rhog = rho;
    *eg   = e;
    *ug   = -u;
    *vg   = v;
    *wg   = w;
    break;
  case 2:	       /* OUTFLOW */
    *rhog = rho;
    *eg   = e;
    *ug   = u;
    *vg   = v;
    *wg   = w;
    break;
    /* problem-specific conditions */
  case 3:
      {
      real xi;
      xi = SIGMASLOPE+1+FLARINGINDEX; //.+FLARINGINDEX;
      *rhog = rho*pow(xg/x,-xi);
      *eg   = e;
      *ug   = 0.0; //colatitude
      *vg   = v; //(v+OMEGAFRAME)*pow(x/xg,1.5)-OMEGAFRAME;//v; //azimuth
      *wg   = 0.0; //radius
      }
      break;
  case 10 : // simple keplerian with dim azim < dim z or dim colat
    *rhog = rho;
    *eg = e;
    *ug = u;
    *vg = (v-OMEGAFRAME)*pow(x/xg,1.5)+OMEGAFRAME;
    *wg = w;
    break;
  case 11 : // keplerian 2-3D : radius (valid for COORDPERMUT = 123 or 213 or 231)
    if (predictive) break;
    *rhog = keplerian_boundary(11,_density_,rho,x,xg,yg,zg);
    *eg = keplerian_boundary(11,_energy_,e,x,xg,yg,zg);
    *ug = keplerian_boundary(11,_vrad_,u,x,xg,yg,zg);
    *vg = keplerian_boundary(11,_vazimuth_,v,x,xg,yg,zg);
    *wg = 0;			/* Colatitude */
    break;
  case 12 : // keplerian 3D : colatitude (valid COORDPERMUT = 213)
    if (predictive) break;
    *rhog = keplerian_boundary(12,_density_,rho,x,xg,yg,zg);
    *eg = keplerian_boundary(12,_energy_,e,x,xg,yg,zg);
    *ug = -u;			/* Velocity in colatitude */
    *wg = w;			/* Velocity in radius */
    *vg = keplerian_boundary(12,_vazimuth_,v,x,xg,yg,zg);
    break;
    /* unknown conditions */
  default:
    prs_error ("Unknow boundary condition code at line %ld in %s",\
	       line, GRIDFILE);
  }
  return Sym;
}

void TrueBC_fp (fluid)
     FluidPatch *fluid;
{
  long gncell[3], stride[3], dim, side, dim1, dim2;
  long i_in, i_gh, j, k, m_in, m_gh, bc, offside;
  real *dens, *vel[3], *energy, *center[3];
  real foo1=0.0, foo2=0.0, foo3=0.0;
  real vp=0.0, *vpg, vt1=0.0, *vtg1, vt2=0.0, *vtg2;
  /* The above initializations are necessary, otherwise when
     'boundary' is first called below to test the symmetry properties,
     a floating point error can be generated (especially on HP/SC) */
  real x, xg;			/* See use of x & xg in the boundary
  function above */
  real yg, zg;
  tGrid_CPU *cpug;
  boolean Sym;
  cpug = fluid->desc;
  dens = fluid->Density->Field;
  energy = fluid->Energy->Field;
  vtg1 = &foo1;
  vtg2 = &foo2;
  for (j = 0; j < 3; j++) {
    vel[j] = fluid->Velocity->Field[j];
    center[j] = fluid->desc->Center[j];
  }
  getgridsize (cpug, gncell, stride);
  for (dim = 0; dim < NDIM; dim++) {
    for (side = INF; side <= SUP; side++) {
      if ((bc = cpug->iface[dim][side]) > 0) {	/* True boundary conditions */
	Sym = boundary (vt1,vt2,vp,vt1,vt2,&foo1,&foo2,&foo1,&foo2, \
			&foo3,1.0,1.0,1.0,1.0,bc,cpug->Parent->linenumber, FALSE);
	/* The above line is just meant to know the symmetry
	   properties of the boundary condition */
	dim1 = (dim == 0);
	dim2 = 2-(dim == 2);
	offside = (side == INF ? 0 : gncell[dim]-Nghost[dim]);
	for (i_gh = offside; i_gh < offside+Nghost[dim]; i_gh++) {
	  i_in = (2*Nghost[dim]-1)-i_gh;
	  if (side == SUP) i_in += 2*(gncell[dim]-2*Nghost[dim]);
	  if (!Sym)
	    i_in = (side == INF ? Nghost[dim] : gncell[dim]-Nghost[dim]-1);
	  for (j = 0; j < gncell[dim1]; j++) {
	    for (k = 0; k < gncell[dim2]; k++) {
	      m_in = i_in*stride[dim]+j*stride[dim1]+k*stride[dim2];
	      m_gh = i_gh*stride[dim]+j*stride[dim1]+k*stride[dim2];
	      vp = vel[dim][m_in];
	      vpg=&vel[dim][m_gh];
	      x = center[dim][m_in];
	      xg= center[dim][m_gh];
	      yg= center[dim1][m_gh];
	      zg= center[dim2][m_gh];
	      if (NDIM > 1) {
		vt1 = vel[dim1][m_in];
		vtg1=&vel[dim1][m_gh];
	      }
	      if (NDIM > 2) {
		vt2 = vel[dim2][m_in];
		vtg2=&vel[dim2][m_gh];
	      }


          boundary (dens[m_in],energy[m_in],vp,vt1,vt2,&dens[m_gh],\
  			&energy[m_gh],vpg,vtg1,vtg2,x,xg,yg,zg, \
  			bc,cpug->Parent->linenumber, FALSE);


	    }
	  }
	}
      }
    }
  }
}

void TrueBC (level)
     long level;
{
  tGrid_CPU *item;
  FluidPatch *fluid;
  item = Grid_CPU_list;
  while (item != NULL) {
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      fluid = item->Fluid;
      while (fluid != NULL) {
      	TrueBC_fp (fluid);
      	fluid = fluid->next;
      }
    }
    item = item->next;
  }
}
