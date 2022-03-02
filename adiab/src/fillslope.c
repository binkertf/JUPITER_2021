#include "jupiter.h"

void FillSlopes () {
  long i, j, k, l, ll, m, n, mp[3], mm[3], ind[3], ip1, ip2, w1, w2, d1, d2;
  FluidWork *fw;
  real *edges[3], *center[3], *rho, *v[3];
  real vp, vm, vc, rhoc, rhop, rhom, *energy;
  real ec, ep, em;
  real *srho[3], *sv[3][3], invdx[3], *invm[3][2], *met[3][2], *se[3];
  real fdx[3];
  real *inter[3];
  long gncell[3], stride[3], ngh[3];
  long gncells[3], strides[3];
  real *rho0=NULL;
  real *ene0=NULL;
  long Size;
  fw = CurrentFluidPatch;
  getgridsizes (fw->desc, gncell, stride, gncells, strides);
  rho = fw->Density;
  Size = gncell[0]*gncell[1]*gncell[2];
  if (HydroStaticEnforce) {
    rho0 = fw->Fluid->Rho_eq_c->field;
    if (!Isothermal)
      ene0 = fw->Fluid->Ene_eq_c->field;
  }
	if (!Isothermal)
  	energy = fw->Energy;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    edges[i] = fw->desc->Edges[i];
    invm[i][0] = fw->desc->InvMetric[i][0];
    invm[i][1] = fw->desc->InvMetric[i][1];
    met[i][0] = fw->desc->Metric[i][0];
    met[i][1] = fw->desc->Metric[i][1];
    center[i] = fw->desc->Center[i];
    srho[i] = fw->SlopeDensity[i];
		if (!Isothermal)
    	se[i]  = fw->SlopeEnergy[i];
    inter[i] = fw->desc->InterSurface[i];
    for (j = 0; j < NDIM; j++) {
      sv[i][j] = fw->SlopeVelocity[i][j];
    }
  }
  for (l = 0; l < NDIM; l++)
	if (!Isothermal){
		for (i = 0; i < Size; i++)
      se[l][i] = 0.0;
	}
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (n = 0; n < NDIM; n++) {
    ip1 = (n == 0);
    ip2 = 2 - (n == 2);
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      ind[2] = k;
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	ind[1] = j;
	for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	  ind[0] = i;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  for (ll = 0; ll < NDIM; ll++) {
	    mp[ll] = m+stride[ll];
	    mm[ll] = m-stride[ll];
	  }
	  invdx[0] = 1./(edges[0][i+1]-edges[0][i]);
	  fdx[0] = invm[0][0][j] * invm[0][1][k];
	  invdx[0] *= fdx[0];
	  if (NDIM > 1) {
	    invdx[1] = 1./(edges[1][j+1]-edges[1][j]);
	    fdx[1] = invm[1][0][i] * invm[1][1][k];
	    invdx[1] *= fdx[1];
	  }
	  if (NDIM > 2) {
	    invdx[2] = 1./(edges[2][k+1]-edges[2][k]);
	    fdx[2] = invm[2][0][i] * invm[2][1][j];
	    invdx[2] *= fdx[2];
	  }
	  for (l = 0; l < NDIM; l++) {
	    if (n == 0) {		/* Slopes of scalar quantities are evaluated only once */
	      rhoc = rho[m];
	      rhop = rho[mp[l]];
	      rhom = rho[mm[l]];
	      srho[l][m] = invdx[l]*__TVDslope(rhop-rhoc,rhoc-rhom);
	      if (!Isothermal) {
		ec = energy[m];
		ep = energy[mp[l]];
		em = energy[mm[l]];
		se[l][m] = invdx[l]*__TVDslope(ep-ec,ec-em);
	      }
	    }
	    w1 = ind[ip1];
	    w2 = ind[ip2];
	    d1 = (ip1 == l);
	    d2 = (ip2 == l);
	    vc = v[n][m]*met[n][0][w1]*met[n][1][w2]; /* We enforce having linear velocities */
	    vp = v[n][mp[l]]*met[n][0][w1+d1]*met[n][1][w2+d2];
	    vm = v[n][mm[l]]*met[n][0][w1-d1]*met[n][1][w2-d2];
	    sv[l][n][m] = invdx[l]*__TVDslope(vp-vc,vc-vm);
	  }
	}
      }
    }
  }
}





void FillSlopes_Dust () {
  long i, j, k, l, ll, m, n, mp[3], mm[3], ind[3], ip1, ip2, w1, w2, d1, d2;
  FluidWork *fw;
  real *edges[3], *center[3], *rho, *v[3];
  real vp, vm, vc, rhoc, rhop, rhom, *energy;
  //  real ec, ep, em;
  real *srho[3], *sv[3][3], invdx[3], *invm[3][2], *met[3][2], *se[3];
  real fdx[3], fvel[3], invf;
  real *inter[3];
  long gncell[3], stride[3], ngh[3];
  long gncells[3], strides[3], ms, msp, msm;
  real *rho0=NULL;
  fw = CurrentFluidPatch;
  getgridsizes (fw->desc, gncell, stride, gncells, strides);
  rho = fw->Density;
  if (HydroStaticEnforce)
    rho0 = fw->Fluid->Rho_eq_c->field;
  energy = fw->Energy;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    edges[i] = fw->desc->Edges[i];
    invm[i][0] = fw->desc->InvMetric[i][0];
    invm[i][1] = fw->desc->InvMetric[i][1];
    met[i][0] = fw->desc->Metric[i][0];
    met[i][1] = fw->desc->Metric[i][1];
    center[i] = fw->desc->Center[i];
    srho[i] = fw->SlopeDensity[i];
    se[i]  = fw->SlopeEnergy[i];
    inter[i] = fw->desc->InterSurface[i];
    for (j = 0; j < NDIM; j++) {
      sv[i][j] = fw->SlopeVelocity[i][j];
    }
  }
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (n = 0; n < NDIM; n++) {
    ip1 = (n == 0);
    ip2 = 2 - (n == 2);
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      ind[2] = k;
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	ind[1] = j;
	for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	  ind[0] = i;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  for (ll = 0; ll < NDIM; ll++) {
	    mp[ll] = m+stride[ll];
	    mm[ll] = m-stride[ll];
	  }
	  invdx[0] = 1./(edges[0][i+1]-edges[0][i]);
	  fdx[0] = invm[0][0][j] * invm[0][1][k];
	  invdx[0] *= fdx[0];
	  if (NDIM > 1) {
	    invdx[1] = 1./(edges[1][j+1]-edges[1][j]);
	    fdx[1] = invm[1][0][i] * invm[1][1][k];
	    invdx[1] *= fdx[1];
	  }
	  if (NDIM > 2) {
	    invdx[2] = 1./(edges[2][k+1]-edges[2][k]);
	    fdx[2] = invm[2][0][i] * invm[2][1][j];
	    invdx[2] *= fdx[2];
	  }
	  for (l = 0; l < NDIM; l++) {
	    if (n == 0) {		/* Slopes of scalar quantities are evaluated only once */
	      rhoc = rho[m];
	      rhop = rho[mp[l]];
	      rhom = rho[mm[l]];
	      if (CorrHydroStat[l]) {
		ms = i*strides[0]+j*strides[1]+k*strides[2];
		msm = ms-strides[l];
		msp = ms+strides[l];
		rhoc = rhoc-rho0[ms];
		rhop = rhop-rho0[msp];
		rhom = rhom-rho0[msm];
	      }
	      srho[l][m] = invdx[l]*__TVDslope(rhop-rhoc,rhoc-rhom);
	      //	  if (!Isothermal) {
	      //	    ec = energy[m];
	      //	    ep = energy[mp[l]];
	      //	    em = energy[mm[l]];
	      //	    se[l][m] = invdx[l]*__TVDslope(ep-ec,ec-em);
	      //	  }
	    }
	    invf = 1.0;
	    w1 = ind[ip1];
	    w2 = ind[ip2];
	    d1 = (ip1 == l);
	    d2 = (ip2 == l);
	    vc = v[n][m];
	    vp = v[n][mp[l]];
	    vm = v[n][mm[l]];
	    fvel[0] = met[n][0][w1] * met[n][1][w2];
	    fvel[1] = met[n][0][w1+d1] * met[n][1][w2+d2];
	    fvel[2] = met[n][0][w1-d1] * met[n][1][w2-d2];
	    if (fvel[0] != fvel[1]) { /* Means that the velocity considered is an angular one */
	      /* in which case we consider rather the associated momentum */
	      /* One test only suffices. This happens only in a transverse case (n != l), */
	      /* which is why  we do not need consider  this in a zone
		 splitting case */
	      fvel[0] *= fvel[0];
	      fvel[1] *= fvel[1];
	      fvel[2] *= fvel[2];
	      invf = fdx[n];
	    }
	    vc *= fvel[0];
	    vp *= fvel[1];
	    vm *= fvel[2];
	    sv[l][n][m] = invf*invdx[l]*__TVDslope(vp-vc,vc-vm);
	  }
	}
      }
    }
  }
}