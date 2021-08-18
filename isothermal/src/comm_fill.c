#include "jupiter.h"

void ExecCommSame (lev)
     long lev;
{
  int fieldtype[20]={_Velocity_, _Density_, _Energy_};
  if (Isothermal)
    ExecCommSameVar (lev, 2, fieldtype);
  else
    ExecCommSameVar (lev, 3, fieldtype);
}

void ExecCommSameVar (long lev, long nb, int *fieldtype)
{
  /* Execute intra-level communications, for level lev. These are
     GHOST type communications */
  Communicator *com;
  long i,j,k,l,m,n,nvar,imin[3],imax[3],stride[3],d;
  int field;
  FluidPatch *fluid;
  real *source[80];
  com = ComListGhost;
  nvar = nbvariables (nb, fieldtype);
  while (com != NULL) {
    if ((com->dest_level == lev) && (com->src_level == lev)) {
      if (com->CPU_src == CPU_Rank) {
	l=0;
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  for (i = 0; i < nb; i++) {
	    field = fieldtype[i];
	    d = dimfield(field);
	    for (j = 0; j < d; j++)
	      source[l++] = fluid->Ptr[field+j];
	  }
	  fluid = fluid->next;
	}
	for (l = 0; l < 3; l++) {
	  imin[l] = com->imin_src[l];
	  imax[l] = com->imax_src[l];
	  stride[l] = com->srcg->stride[l];
	}
	n=0;
	for (l=0; l < nvar*NbFluids; l++) {
	  for (k = imin[2]; k < imax[2]; k++) {
	    for (j = imin[1]; j < imax[1]; j++) {
	      for (i = imin[0]; i < imax[0]; i++) {
		m = i*stride[0]+j*stride[1]+k*stride[2];
		com->buffer[n++] = source[l][m];
	      }
	    }
	  }
	}
	if (com->size*nvar*NbFluids != n)
	  prs_error ("Internal error: communicator size mismatch error (lev2lev)");
      }
    }
    com = com->next;
  }
  ExecComm (lev,lev,GHOST,nvar,nb,fieldtype);
}

void ExecCommUp (lev)
     long lev;
{
  int fieldtype[20]={_Density_, _Velocity_, _Energy_};
  if (Isothermal)
    __ExecCommUpVar (lev, 2, fieldtype);
  else
    __ExecCommUpVar (lev, 3, fieldtype);
}

void ExecCommUpVar (long lev, long nb, int *fieldtype)
{
  /* Execute communications from level lev to finer level lev+1. These
     are GHOST type communications. In this first version we use a direct injection */
  Communicator *com;
  long i,j,k,l,ms,n,nvar,imind[3],imaxd[3],d;
  long imins[3], strides[3],is,js,ks,idx[3];
  long comp[80], sqz[]={0,0,0};
  real *source[80], *radius, *colatitude, s;
  real rad, col;
  int field;
  FluidPatch *fluid;
  com = ComListGhost;
  nvar = nbvariables (nb, fieldtype);
  while (com != NULL) {
    if ((com->dest_level == lev+1) && (com->src_level == lev) && (com->type == GHOST)) {
      if (com->CPU_src == CPU_Rank) {
	l=0;
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  for (i = 0; i < nb; i++) {
	    field = fieldtype[i];
	    d = 1;
	    if (field == _Velocity_)
	      d = NDIM;
	    for (j = 0; j < d; j++) {
	      comp[l] = _other_;
	      if (field == _Density_) comp[l] = _density_;
	      if (field == _Energy_) comp[l] = _energy_;
	      if ((field == _Velocity_) && (j == _AZIM_)) comp[l] = _vazimuth_;
	      source[l++] = fluid->Ptr[field+j];
	    }
	  }
	  fluid = fluid->next;
	}
	radius     = com->srcg->Fluid->desc->Center[_RAD_];
	colatitude = com->srcg->Fluid->desc->Center[_COLAT_];
	for (l = 0; l < 3; l++) {
	  imind[l] = com->imin_dest[l];
	  imaxd[l] = com->imax_dest[l];
	  imins[l] = com->imin_src[l];
	  strides[l] = com->srcg->stride[l];
	  if (HydroStaticReady)
	    sqz[l] = com->srcg->Fluid->Rho_eq_c->stride[l];
	}
	n=0;
	for (l=0; l<nvar*NbFluids; l++) {
	  for (k = imind[2]; k < imaxd[2]; k++) {
	    for (j = imind[1]; j < imaxd[1]; j++) {
	      for (i = imind[0]; i < imaxd[0]; i++) {
		is = (i-imind[0])/(Refine[0] ? 2:1)+imins[0];
		js = (j-imind[1])/(Refine[1] ? 2:1)+imins[1];
		ks = (k-imind[2])/(Refine[2] ? 2:1)+imins[2];
		idx[2] = ks;
		idx[1] = js;
		idx[0] = is;
		ms = is*strides[0]+js*strides[1]+ks*strides[2];
		s = source[l][ms];
		rad = radius[ms];
		col = colatitude[ms];
		s = comm_adapt (s, comp[l], ms, idx, sqz, com->srcg, +1, l/nvar);
		com->buffer[n++] = s;
	      }
	    }
	  }
	}
	if (com->size*nvar*NbFluids != n)
	  prs_error ("Internal error: communicator size mismatch error");
      }
    }
    com = com->next;
  }
  ExecComm (lev,lev+1,GHOST,nvar,nb,fieldtype);
}

void ExecCommUpVarLIL (long lev, long nb, int *fieldtype) /* With slope limiter */
{
  /* Execute communications from level lev to finer level lev+1. These
     are GHOST type communications. In this version we use a (multi)-linear interpolation */
  Communicator *com;
  long g,h,i,j,k,l,ms,mc,n,nvar,imind[3],imaxd[3],stridec[3], imaxs[3];
  long imins[3], strides[3], is[3], iid[3], size, sizec[3], msm, msp, d, isl[3];
  real *source[80], s, src, srcp, srcm, loc_slope[3], ic[3];
  int field;
  long comp[20], sqz[]={0,0,0};
  tGrid_CPU *desc;
  FluidPatch *fluid;
  com = ComListGhost;
  nvar = nbvariables (nb, fieldtype);
  while (com != NULL) {
    if ((com->dest_level == lev+1) &&\
	(com->src_level == lev) &&\
	(com->type == GHOST)) {
      if (com->CPU_src == CPU_Rank) {
	l=0;
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  for (i = 0; i < nb; i++) {
	    field = fieldtype[i];
	    d = 1;
	    if (field == _Velocity_)
	      d = NDIM;
	    for (j = 0; j < d; j++) {
	      comp[l] = _other_;
	      if (field == _Density_) comp[l] = _density_;
	      if (field == _Energy_) comp[l] = _energy_;
	      if ((field == _Velocity_) && (j == _AZIM_))
		comp[l] = _vazimuth_;
	      source[l++] = fluid->Ptr[field+j];
	    }
	  }
	  fluid = fluid->next;
	}
	size = com->size;
	for (l = 0; l < 3; l++) {
	  imind[l] = com->imin_dest[l];
	  imaxd[l] = com->imax_dest[l];
	  imins[l] = com->imin_src[l];
	  imaxs[l] = com->imax_src[l];
	  strides[l] = com->srcg->stride[l];
	  if (HydroStaticReady)
	    sqz[l] = com->srcg->Fluid->Rho_eq_c->stride[l];
	  sizec[l] = imaxd[l]-imind[l]; /* Communicator size */
	  stridec[l] = (l > 0 ? stridec[l-1]*sizec[l-1] : 1); /* Communicator stride */
	}
	n=0;
	desc = com->srcg;
	for (l=0; l<nvar*NbFluids; l++) {
	  for (k = imins[2]; k < imaxs[2]; k++) {
	    is[2] = k;
	    for (j = imins[1]; j < imaxs[1]; j++) {
	      is[1] = j;
	      for (i = imins[0]; i < imaxs[0]; i++) {
		is[0] = i;
		ms = 0;
		for (h = 0; h < 3; h++)
		  ms += is[h]*strides[h];
		src = source[l][ms]; /* Source value at coarse point */
		src = comm_adapt (src, comp[l], ms, is, sqz, desc, +1, l/nvar);
		/* We evaluate the three limited slopes at current source (coarse) point */
		/* The following loop is correct even if NDIM<3*/
		for (h = 0; h < 3; h++) {
		  if (Refine[h]) {
		    msp = ms+strides[h];
		    srcp = source[l][msp];
		    for (g = 0; g < 3; g++)
		      isl[g] = (g == h ? is[g]+1 : is[g]);
		    srcp = comm_adapt (srcp, comp[l], msp, isl, sqz, desc, +1, l/nvar);
		      msm = ms-strides[h];
		    srcm = source[l][msm];
		    for (g = 0; g < 3; g++)
		      isl[g] = (g == h ? is[g]-1 : is[g]);
		    srcm = comm_adapt (srcm, comp[l], msm, isl, sqz, desc, +1, l/nvar);
		      loc_slope[h] = .25 * TVDslope (srcp-src,src-srcm);
		  } else {
		    loc_slope[h] = 0.0;
		  }
		}
		for (iid[2] = 0; iid[2] <= Refine[2]; iid[2]++) { /* We scan the 1 to 8 zones */
		  for (iid[1] = 0; iid[1] <= Refine[1]; iid[1]++) { /* covered by the coarse zone */
		    for (iid[0] = 0; iid[0] <= Refine[0]; iid[0]++) {
		      /* We evaluate the target (dest) coordinate within communicator */
		      for (h = 0; h < 3; h++)
			ic[h] = (is[h]-imins[h])*(Refine[h] ? 2 : 1)+iid[h];
		      /* ic stands for index in communicator */
		      mc = 0;
		      for (h = 0; h < 3; h++)
			mc += ic[h]*stridec[h];
		      s = src;
		      for (h = 0; h < 3; h++)  {
			if (Refine[h])
			  s += loc_slope[h]*(iid[h]*2-1);
		      }
		      com->buffer[mc+size*l] = s;
		      n++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	if (size*nvar*NbFluids != n)
	  prs_error ("Internal error: communicator size mismatch error");
      }
    }
    com = com->next;

  }
  ExecComm (lev,lev+1,GHOST,nvar,nb,fieldtype);
}

void ExecCommDownMean (lev)
     long lev;
{
  /* Execute communications from level lev to coarser level lev-1, of
     type MEAN */
  Communicator *com;
  long i,j,k,l,ms,n,nvar,imind[3],imaxd[3], folded, ii[3], sqz[3];
  long imins[3], imaxs[3], stridec[3], strides[3],ic,jc,kc;
  real *source[80], src;
  tGrid_CPU *desc;
  long d,nb=2,comp[20];
  int field, fieldtype[20]={_Density_, _Velocity_, _Energy_};
  FluidPatch *fluid;
  if (!Isothermal) nb++;
  com = ComListMean;
  nvar = 1+NDIM;		/* Density (or mass flux) + velocities (or momentum flux) */
  if (!Isothermal) nvar++;	/* Energy or energy flux */
  while (com != NULL) {
    if ((com->dest_level == lev-1) && (com->src_level == lev) && (com->type == MEAN)) {
      if (com->CPU_src == CPU_Rank) {
	l=0;
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  for (i = 0; i < nb; i++) {
	    field = fieldtype[i];
	    d = 1;
	    if (field == _Velocity_)
	      d = NDIM;
	    for (j = 0; j < d; j++) {
	      comp[l] = _other_;
	      if (field == _Density_) comp[l] = _density_;
	      if (field == _Energy_) comp[l] = _energy_;
	      if ((field == _Velocity_) && (j == _AZIM_)) comp[l] = _vazimuth_;
	      source[l++] = fluid->Ptr[field+j];
	    }
	  }
	  fluid = fluid->next;
	}
	desc = com->srcg;
	for (l = 0; l < 3; l++) {
	  imind[l] = com->imin_dest[l];
	  imaxd[l] = com->imax_dest[l];
	  imins[l] = com->imin_src[l];
	  imaxs[l] = com->imax_src[l];
	  stridec[l] = (l > 0 ? stridec[l-1]*\
			(imaxd[l-1]-imind[l-1]) : 1);
	  strides[l] = desc->stride[l];
	  if (HydroStaticReady)
	    sqz[l] = desc->Fluid->Rho_eq_c->stride[l];
	}
	/* We first reset the buffer */
	for (n=0; n<com->size*nvar*NbFluids; n++)
	  com->buffer[n] = 0.0;
	folded = (Refine[0] ? 2:1) * (Refine[1] ? 2:1) * (Refine[2] ? 2:1);
	for (l=0; l<nvar*NbFluids; l++) {
	  for (k = imins[2]; k < imaxs[2]; k++) {
	    for (j = imins[1]; j < imaxs[1]; j++) {
	      for (i = imins[0]; i < imaxs[0]; i++) {
		ms = i*strides[0]+j*strides[1]+k*strides[2];
		ii[0] = i;
		ii[1] = j;
		ii[2] = k;
		ic = (i-imins[0])/(Refine[0] ? 2:1);
		jc = (j-imins[1])/(Refine[1] ? 2:1);
		kc = (k-imins[2])/(Refine[2] ? 2:1);
		n = ic*stridec[0]+jc*stridec[1]+kc*stridec[2];
		src = source[l][ms];
		src = comm_adapt (src, comp[l], ms, ii, sqz, desc, +1, l/nvar);
		com->buffer[l*com->size+n] += src/(real)folded;
	      }
	    }
	  }
	}
      }
    }
    com = com->next;
  }
  fieldtype[0] = _Density_;
  fieldtype[1] = _Velocity_;
  fieldtype[2] = _Energy_;	/* not useful */
  ExecComm (lev,lev-1,MEAN,nvar,2,fieldtype);
}

void ExecCommDownFlux (lev)
     long lev;
{
  /* Execute communications from level lev to coarser level lev-1, of
     type FLUX */
  Communicator *com;
  long nvar, dim, side, i, imin[3], imax[3], dr[3];
  long size[3], csized[3], csizes[3];
  long dimp1, dimp2, n, l, j, m, buf_size;
  real *source[80];
  FluidPatch *fluid;
  com = ComListFlux;
  nvar = 2+NDIM;		/* Density (or mass flux) + velocities (or momentum flux) */
				/* and interface pressure at the outer faces */
  if (!Isothermal) nvar++;	/* Energy of energy flux */
  while (com != NULL) {
    if ((com->dest_level == lev-1) && (com->src_level == lev) && (com->type == FLUX)) {
      if (com->CPU_src == CPU_Rank) {
	dim  = com->facedim;
	side = com->faceside;
	for (i = 0; i < 3; i++) {
	  dr[i] = (Refine[i] ? 2:1);
	  imin[i] = com->imin_src[i];
	  imax[i] = com->imax_src[i];
	  size[i] = com->srcg->ncell[i];
	  csizes[i] = com->imax_src[i]-com->imin_src[i];
	  csized[i] = com->imax_dest[i]-com->imin_dest[i];
	  if ((csizes[i]/dr[i] != csized[i]) && (i != dim))
	    prs_error ("Flux communicator size internal error.");
	}
	l = 0;
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  source[l++] = fluid->MassFlux->Flux[dim][side];
    //source[l++] = fluid->DiffFlux->Flux[dim][side];
	  for (i = 1; i <= NDIM; i++)
	    source[l++] = fluid->MomentumFlux[i-1]->Flux[dim][side];
	  source[l++] = fluid->Pressure->Pressure[dim][side];
	  if (!Isothermal)
	    source[l++] = fluid->EnergyFlux->Flux[dim][side];
	  fluid = fluid->next;
	}
	dimp1 = (dim == 0);
	dimp2 = 2 - (dim == 2);
	buf_size = csized[dimp1]*csized[dimp2];
	for (i = 0; i < buf_size*nvar*NbFluids; i++)
	  com->buffer[i] = 0.0;
	for (l = 0; l < nvar*NbFluids; l++) {
	  for (i = imin[dimp1]; i < imax[dimp1]; i++) {
	    for (j = imin[dimp2]; j < imax[dimp2]; j++) {
	      m = (i-Nghost[dimp1])+(j-Nghost[dimp2])*size[dimp1];
	      n = (i-imin[dimp1])/dr[dimp1]+\
		(j-imin[dimp2])/dr[dimp2]*csized[dimp1];
	      com->buffer[n+l*buf_size] += source[l][m];
	    }
	  }
	}
      }
    }
    com = com->next;
  }
  ExecCommFlux (lev);
}
