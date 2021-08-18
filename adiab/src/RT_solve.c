#include "jupiter.h"

real FLD (real erre) {
  real lambda;
  if (erre <=2.) {
    lambda = 2./(3.+sqrt(9. +10.*erre*erre));
  } else {
    lambda = 10./(10.*erre+9.+sqrt(81.+180.*erre));
  }
  return lambda;
}



void ComputeDiffusionCoefficients () {
  FluidWork *fw;
  real *radint, *radius, *colint, *azimint, *colat;
  real *dens, *opar, *diffcoeff, *enrad;
  real eip, eim, ejp, ejm, ehp, ehm, erre, coeff, dcolat;
  long gncell[3], stride[3];
  long i,j,h,l,lip,lim,ljp,ljm,lhp,lhm;
  real steprad, stepazim, stepcol, rmed, invdiffr, invdcolat, invdazim;
  real graderrad, gradercol, graderazi;
  real coefr, coeft, coefp;
  fw = CurrentFluidPatch;
  diffcoeff = fw->Diffcoeff;
  dens      = fw->Density;
  opar      = fw->OpaR;
  enrad     = fw->EnergyRad;
  getgridsize (fw->desc, gncell, stride);
  
  radint = fw->desc->Edges[_RAD_];
  radius = fw->desc->Center[_RAD_];
  colat  = fw->desc->Center[_COLAT_];
  colint = fw->desc->Edges[_COLAT_];
  azimint = fw->desc->Edges[_AZIM_];

  steprad = radint[1]-radint[0];
  stepcol = colint[1]-colint[0];
  stepazim = azimint[1]-azimint[0];


// Loop on every cells to compute diffusion coefficients:
  for(h = 0; h < gncell[2]; h++) {
    for (i = 0; i < gncell[1]; i++) {
      rmed = .5*(radint[i]+radint[i+1]);
      invdiffr = 1./(radint[i+1]-radint[i]);
      for (j = 0; j < gncell[0]; j++) {
	l = j+i*stride[1]+h*stride[2];
	dcolat = stepcol*radius[l];
	invdcolat = 1.0/dcolat;
	invdazim = 1./(rmed*sin(colat[l])*stepazim);
	lim = l-stride[1];
	ljm = l-1;
	lhm = l-stride[2];
	lip = l+stride[1];
	ljp = l+1;
	lhp = l+stride[2];
	coeff = 1.0/(opar[l]*dens[l]);
	eip   = (i < gncell[1]-1 ? enrad[lip] : enrad[l]);
	coefr = (i < gncell[1]-1 ? .5         : 1.0     );
	eim   = (i > 0           ? enrad[lim] : enrad[l]);
	coefr = (i > 0           ? .5         : 1.0     );
	ejp   = (j < gncell[0]-1 ? enrad[ljp] : enrad[l]);
	coeft = (j < gncell[0]-1 ? .5         : 1.0     );
	ejm   = (j > 0           ? enrad[ljm] : enrad[l]);
	coeft = (j > 0           ? .5         : 1.0     );
	ehp   = (h < gncell[2]-1 ? enrad[lhp] : enrad[l]);
	coefp = (h < gncell[2]-1 ? .5         : 1.0     );
	ehm   = (h > 0           ? enrad[lhm] : enrad[l]);
	coefp = (h > 0           ? .5         : 1.0     );
	graderrad = coefr*(eip-eim)*invdiffr;
	graderazi = coeft*(ejp-ejm)*invdazim;
	gradercol = coefp*(ehp-ehm)*invdcolat; 
	erre = coeff*sqrt(graderrad*graderrad+graderazi*graderazi+gradercol*gradercol)/enrad[l];
	diffcoeff[l] = FLD(erre)*coeff*C0;
      }
    }
  }

  ExecCommSameOneField (fw->desc->level, diffcoeff); 
}

void ComputeMatrixElements(real dt)
{
  FluidWork *fw;
  long size[3], gncell[3], stride[3];
  long i, j,h,k, l,lip,lim,ljp,ljm,lhp,lhm, nr, ns, ni;
  real *dens,*opar, *opap,*eta1,*eta2,*temper,*enrad,*diffcoeff;
  real *colint, *radint, *colat;
  real *aar,*aap,*aat,*ccr,*cct,*ccp,*bb,*rhs;
  real temper3,sinphip,sinphim,cosphip,cosphim;
  real rsup,rinf,rmed,invdiffr,dxcol,invdxcol,stepcol,stepazim,invdxazi;
  real denom;

  denom = .5*dt;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  
  for (i = 0; i < 3; i++)
    size[i] = fw->desc->gncell[i];
  
  nr = size[1]; // if coordpermut=213 (first coord is azimuth, 2nd rad, 3rd colat)
  ns = size[0]; // azimuth
  ni = size[2]; // co-lat
  
  dens   =  fw->Density; 
  opar   =  fw->OpaR;
  opap   =  fw->OpaP;
  temper = fw->Temperature;
  enrad  = fw->EnergyRad;
  
  eta1 = fw->Eta1;
  eta2 = fw->Eta2;
  diffcoeff = fw->Diffcoeff;
  aar = fw->Aarmatrix;
  aat = fw->Aatmatrix;
  aap = fw->Aapmatrix;
  ccr = fw->Ccrmatrix;
  cct = fw->Cctmatrix;
  ccp = fw->Ccpmatrix;
  bb  = fw->Bbmatrix;
  rhs = fw->Rhsmatrix;
  
  radint  = fw->desc->Edges[_RAD_];
  colint  = fw->desc->Edges[_COLAT_];
  colat   = fw->desc->Center[_COLAT_];
  stepcol = colint[1]-colint[0];
  stepazim= fw->desc->Edges[_AZIM_][1]-fw->desc->Edges[_AZIM_][0];
  
  for (h = 1; h < gncell[2]-1; h++) {
    sinphip = sin(colint[h+1]);
    sinphim = sin(colint[h]);
    cosphip = cos(colint[h+1]);
    cosphim = cos(colint[h]);
    for (i = 1; i < gncell[1]-1; i++) {
      k = h*stride[2]+i*stride[1]; 
      rsup = radint[i+1];
      rinf = radint[i];
      rmed = .5*(rsup+rinf);
      invdiffr=1./(rsup-rinf);
      for (j = 1; j < gncell[0]-1; j++) {
	l = j+i*stride[1]+h*stride[2];
        dxcol = stepcol*rmed;
        invdxcol = 1.0/dxcol;
	invdxazi = 1./(rmed*sin(colat[l])*stepazim);
        lim = l-stride[1];
	lip = l+stride[1];
	ljp = l+1;
	ljm = l-1;
	lhp = l+stride[2];
	lhm = l-stride[2];
	ccr[l]  = -denom*(diffcoeff[lip]+diffcoeff[l])*pow(rsup/rmed*invdiffr,2.0);
       	aar[l]  = -denom*(diffcoeff[l]+diffcoeff[lim])*pow(rinf/rmed*invdiffr,2.0);
	cct[l]  = -denom*(diffcoeff[l]+diffcoeff[ljp])*invdxazi*invdxazi;
	aat[l]  = -denom*(diffcoeff[l]+diffcoeff[ljm])*invdxazi*invdxazi;
	ccp[l]  = -denom*(diffcoeff[lhp]+diffcoeff[l])*invdxcol*sinphip/(cosphim-cosphip)/rmed;
        aap[l]  = -denom*(diffcoeff[l]+diffcoeff[lhm])*invdxcol*sinphim/(cosphim-cosphip)/rmed;
	bb[l]   = -aar[l]-ccr[l]-aat[l]-cct[l]-aap[l]-ccp[l];
	temper3 = pow(temper[l],3.0);
	bb[l]   = 1.0+dt*dens[l]*opap[l]*(C0-16.0*SIGMARAD*temper3*eta2[l])+bb[l];
	rhs[l]  = enrad[l]+4.0*dt*dens[l]*opap[l]*SIGMARAD*temper3*(4.0*eta1[l]-3.0*temper[l]);
      }
    }
  }
}
