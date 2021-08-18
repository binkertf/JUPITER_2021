#include "jupiter.h"

void FillBeam  (dim, ii, jj, beam)
     long dim, ii, jj;
     Beam *beam;
{
  long k, gncell[3], stride[3], ip[2], offset, offsets;
  long strides[3], gncells[3];
  real metpar, *met[3][2], *metperp[2];
  long idim, ldim[3], str, len, strs;
  FluidWork *fw;
  fw = CurrentFluidPatch;
  for (k = 0; k < 3; k++) {
    met[k][0] = fw->desc->Metric[k][0];
    met[k][1] = fw->desc->Metric[k][1];
  }
  metperp[0] = beam->metperp[0];
  metperp[1] = beam->metperp[1];
  ldim[0] = dim;		/* ldim : {#dim beam axis, #perp dim 1 beam, #perp dim 2 beam} */
  ldim[1] = ip[0] = (dim == 0);
  ldim[2] = ip[1] = 2-(dim == 2);
  getgridsizes (fw->desc, gncell, stride, gncells, strides);
  offset = ii*stride[ip[0]]+jj*stride[ip[1]];
  offsets= ii*strides[ip[0]]+jj*strides[ip[1]];
  str = stride[dim];
  strs = strides[dim];
  len = gncell[dim];
  metpar = met[dim][0][ii]*met[dim][1][jj];
  if (NDIM > 1) {
    memcpy (metperp[0], met[ip[0]][(ip[1] < dim)], (size_t)(len*sizeof(real)));
    multarray (metperp[0], met[ip[0]][1-(ip[1] < dim)][jj], len);
  }
  if (NDIM > 2) {
    memcpy (metperp[1], met[ip[1]][(ip[0] < dim)], (size_t)(len*sizeof(real)));
    multarray (metperp[1], met[ip[1]][1-(ip[0] < dim)][ii], len);
  }
  if (CorrHydroStat[dim]) {
    memcpystride (beam->HS_cent_rho, fw->Fluid->Rho_eq_c->field+offsets, len, strs);
    memcpystride (beam->HS_int_rho, fw->Fluid->Rho_eq_i[dim]->field+offsets, len, strs);
    memcpystride (beam->cs2i, fw->Fluid->Cs2_i[dim]->field+offsets, len, strs);
  }
  memcpystride (beam->rho, fw->Density+offset, len, str);
  memcpystride (beam->rho_pred, fw->Density_Pred+offset, len, str);
  memcpystride (beam->srcrho, fw->SourceRhoPred+offset, len, str);
  memcpystride (beam->source, fw->SourceVelocity[dim]+offset, len, str);
  memcpystride (beam->rawcoord, fw->desc->Center[dim]+offset, len, str);
  beam->rawcoord1 = *(fw->desc->Center[ldim[1]]+offset);
  beam->radius = *(fw->desc->Center[_RAD_]+offset);
  if (__SPHERICAL && (dim == _AZIM_))
    beam->radius *= sin(*(fw->desc->Center[_COLAT_]+offset));
  if (KEPLERIAN && (dim == _AZIM_)) {
    beam->masscorr1 = fw->Fluid->MassFluxCorrection1[ii+jj*gncell[ldim[1]]];
    beam->masscorr2 = fw->Fluid->MassFluxCorrection2[ii+jj*gncell[ldim[1]]];
    beam->momcorr1 = fw->Fluid->MomentumFluxCorrection1[ii+jj*gncell[ldim[1]]];
    beam->momcorr2 = fw->Fluid->MomentumFluxCorrection2[ii+jj*gncell[ldim[1]]];
    if (!Isothermal) {
      beam->enercorr1 = fw->Fluid->EnergyFluxCorrection1[ii+jj*gncell[ldim[1]]];
      beam->enercorr2 = fw->Fluid->EnergyFluxCorrection2[ii+jj*gncell[ldim[1]]];
    }
  }
  beam->rawcoord2 = *(fw->desc->Center[ldim[2]]+offset);
  memcpystride (beam->u, fw->Velocity[dim]+offset, len, str);
  memcpystride (beam->u_pred, fw->Velocity_Pred[dim]+offset, len, str);
  multarray (beam->u, metpar, len);
  multarray (beam->u_pred, metpar, len);
  memcpystride (beam->center, fw->desc->Center[dim]+offset, len, str);
  multarray (beam->center, metpar, len);
  memcpystride (beam->edge, fw->desc->Edges[dim], len, 1);
  multarray (beam->edge, metpar, len);
  memcpystride (beam->intersurface, fw->desc->InterSurface[dim]+offset, len, str);
  for (idim = 0; idim < NDIM-1; idim++) {
    memcpystride (beam->source_perp[idim], fw->SourceVelocity[ip[idim]]+offset, len, str);
    memcpystride (beam->v_perp[idim], fw->Velocity[ip[idim]]+offset, len, str);
    memcpystride (beam->v_perp_pred[idim], fw->Velocity_Pred[ip[idim]]+offset, len, str);
    arraymult (beam->v_perp[idim], metperp[idim], len);
    arraymult (beam->v_perp_pred[idim], metperp[idim], len);
    beam->MomCorr[idim] = MomentumCorrection[dim][ip[idim]];
    if (beam->MomCorr[idim]) {
      memcpy (beam->MomCorrRatioL[idim],fw->desc->MomCorrRatioL[dim][ip[idim]], (size_t)(len*sizeof(real)));
      memcpy (beam->MomCorrRatioR[idim],fw->desc->MomCorrRatioR[dim][ip[idim]], (size_t)(len*sizeof(real)));
    }
  }
  for (idim = 0; idim < NDIM; idim++) {
    memcpystride (beam->slope_rho[idim], fw->SlopeDensity[ldim[idim]]+offset, len, str);
    if (!Isothermal)
      memcpystride (beam->slope_energy[idim], fw->SlopeEnergy[ldim[idim]]+offset, len, str);
    memcpystride (beam->slope_u[idim], fw->SlopeVelocity[ldim[idim]][ldim[0]]+offset, len, str);
    if (NDIM > 1)
      memcpystride (beam->slope_v_perp[0][idim], fw->SlopeVelocity[ldim[idim]][ldim[1]]+offset, len, str);
    if (NDIM > 2)
      memcpystride (beam->slope_v_perp[1][idim], fw->SlopeVelocity[ldim[idim]][ldim[2]]+offset, len, str);
  }
  for (k = 0; k < len; k++) {
    if (Isothermal)
      beam->cs[k] = sqrt(fw->Energy[offset]);
    else
      beam->cs[k] = fw->Energy[offset];	/* The (unappropriately named)
					   cs field then contains the
					   volumic internal energy */
    offset += str;
  }
  beam->dim[0] = dim;
  beam->dim[1] = ip[0];
  beam->dim[2] = ip[1];
  beam->length = gncell[dim];
  beam->desc = fw->desc;
  beam->true_bc[INF] = fw->desc->iface[dim][INF];
  beam->true_bc[SUP] = fw->desc->iface[dim][SUP];

/*
    if ((dim==1)&&(fw->next!=NULL)){
    if (fw->desc->iface[dim][INF]>0){
      beam->true_bc[INF]=2;
    }
    if (fw->desc->iface[dim][SUP]>0){
      beam->true_bc[SUP]=2;
    }
  }
*/





}


void FillBeam2  (dim, ii, jj, beam)
     long dim, ii, jj;
     Beam *beam;
{
  long k, gncell[3], stride[3], ip[2], offset, offsets;
  long strides[3], gncells[3];
  real metpar, *met[3][2], *metperp[2];
  long idim, ldim[3], str, len, strs;
  FluidWork *fw;
  fw = SecondaryFluidPatch;
  for (k = 0; k < 3; k++) {
    met[k][0] = fw->desc->Metric[k][0];
    met[k][1] = fw->desc->Metric[k][1];
  }
  metperp[0] = beam->metperp[0];
  metperp[1] = beam->metperp[1];
  ldim[0] = dim;		/* ldim : {#dim beam axis, #perp dim 1 beam, #perp dim 2 beam} */
  ldim[1] = ip[0] = (dim == 0);
  ldim[2] = ip[1] = 2-(dim == 2);
  getgridsizes (fw->desc, gncell, stride, gncells, strides);
  offset = ii*stride[ip[0]]+jj*stride[ip[1]];
  offsets= ii*strides[ip[0]]+jj*strides[ip[1]];
  str = stride[dim];
  strs = strides[dim];
  len = gncell[dim];
  metpar = met[dim][0][ii]*met[dim][1][jj];
  if (NDIM > 1) {
    memcpy (metperp[0], met[ip[0]][(ip[1] < dim)], (size_t)(len*sizeof(real)));
    multarray (metperp[0], met[ip[0]][1-(ip[1] < dim)][jj], len);
  }
  if (NDIM > 2) {
    memcpy (metperp[1], met[ip[1]][(ip[0] < dim)], (size_t)(len*sizeof(real)));
    multarray (metperp[1], met[ip[1]][1-(ip[0] < dim)][ii], len);
  }
  if (CorrHydroStat[dim]) {
    memcpystride (beam->HS_cent_rho, fw->Fluid->Rho_eq_c->field+offsets, len, strs);
    memcpystride (beam->HS_int_rho, fw->Fluid->Rho_eq_i[dim]->field+offsets, len, strs);
    memcpystride (beam->cs2i, fw->Fluid->Cs2_i[dim]->field+offsets, len, strs);
  }
  memcpystride (beam->rho, fw->Density+offset, len, str);
  memcpystride (beam->rho_pred, fw->Density_Pred+offset, len, str);
  memcpystride (beam->srcrho, fw->SourceRhoPred+offset, len, str);
  memcpystride (beam->source, fw->SourceVelocity[dim]+offset, len, str);
  memcpystride (beam->rawcoord, fw->desc->Center[dim]+offset, len, str);
  beam->rawcoord1 = *(fw->desc->Center[ldim[1]]+offset);
  beam->radius = *(fw->desc->Center[_RAD_]+offset);
  if (__SPHERICAL && (dim == _AZIM_))
    beam->radius *= sin(*(fw->desc->Center[_COLAT_]+offset));
  if (KEPLERIAN && (dim == _AZIM_)) {
    beam->masscorr1 = fw->Fluid->MassFluxCorrection1[ii+jj*gncell[ldim[1]]];
    beam->masscorr2 = fw->Fluid->MassFluxCorrection2[ii+jj*gncell[ldim[1]]];
    beam->momcorr1 = fw->Fluid->MomentumFluxCorrection1[ii+jj*gncell[ldim[1]]];
    beam->momcorr2 = fw->Fluid->MomentumFluxCorrection2[ii+jj*gncell[ldim[1]]];
    if (!Isothermal) {
      beam->enercorr1 = fw->Fluid->EnergyFluxCorrection1[ii+jj*gncell[ldim[1]]];
      beam->enercorr2 = fw->Fluid->EnergyFluxCorrection2[ii+jj*gncell[ldim[1]]];
    }
  }

  beam->rawcoord2 = *(fw->desc->Center[ldim[2]]+offset);
  memcpystride (beam->u, fw->Velocity[dim]+offset, len, str);
  memcpystride (beam->u_pred, fw->Velocity_Pred[dim]+offset, len, str);
  multarray (beam->u, metpar, len);
  multarray (beam->u_pred, metpar, len);
  memcpystride (beam->center, fw->desc->Center[dim]+offset, len, str);
  multarray (beam->center, metpar, len);
  memcpystride (beam->edge, fw->desc->Edges[dim], len, 1);
  multarray (beam->edge, metpar, len);
  memcpystride (beam->intersurface, fw->desc->InterSurface[dim]+offset, len, str);
  for (idim = 0; idim < NDIM-1; idim++) {
    memcpystride (beam->source_perp[idim], fw->SourceVelocity[ip[idim]]+offset, len, str);
    memcpystride (beam->v_perp[idim], fw->Velocity[ip[idim]]+offset, len, str);
    memcpystride (beam->v_perp_pred[idim], fw->Velocity_Pred[ip[idim]]+offset, len, str);
    arraymult (beam->v_perp[idim], metperp[idim], len);
    arraymult (beam->v_perp_pred[idim], metperp[idim], len);
    beam->MomCorr[idim] = MomentumCorrection[dim][ip[idim]];
    if (beam->MomCorr[idim]) {
      memcpy (beam->MomCorrRatioL[idim],fw->desc->MomCorrRatioL[dim][ip[idim]], (size_t)(len*sizeof(real)));
      memcpy (beam->MomCorrRatioR[idim],fw->desc->MomCorrRatioR[dim][ip[idim]], (size_t)(len*sizeof(real)));
    }
  }
  for (idim = 0; idim < NDIM; idim++) {
    memcpystride (beam->slope_rho[idim], fw->SlopeDensity[ldim[idim]]+offset, len, str);
    if (!Isothermal)
      memcpystride (beam->slope_energy[idim], fw->SlopeEnergy[ldim[idim]]+offset, len, str);
    memcpystride (beam->slope_u[idim], fw->SlopeVelocity[ldim[idim]][ldim[0]]+offset, len, str);
    if (NDIM > 1)
      memcpystride (beam->slope_v_perp[0][idim], fw->SlopeVelocity[ldim[idim]][ldim[1]]+offset, len, str);
    if (NDIM > 2)
      memcpystride (beam->slope_v_perp[1][idim], fw->SlopeVelocity[ldim[idim]][ldim[2]]+offset, len, str);
  }
  for (k = 0; k < len; k++) {
    if (Isothermal)
      beam->cs[k] = sqrt(fw->Energy[offset]);
    else
      beam->cs[k] = fw->Energy[offset];	/* The (unappropriately named)
					   cs field then contains the
					   volumic internal energy */
    offset += str;
  }
  beam->dim[0] = dim;
  beam->dim[1] = ip[0];
  beam->dim[2] = ip[1];
  beam->length = gncell[dim];
  beam->desc = fw->desc;
  beam->true_bc[INF] = fw->desc->iface[dim][INF];
  beam->true_bc[SUP] = fw->desc->iface[dim][SUP];

/*
    if ((dim==1)&&(fw->next!=NULL)){
    if (fw->desc->iface[dim][INF]>0){
      beam->true_bc[INF]=2;
    }
    if (fw->desc->iface[dim][SUP]>0){
      beam->true_bc[SUP]=2;
    }
  }
*/





}


void FillFluxes  (dim, ii, jj, beam)
     long dim, ii, jj;
     Beam *beam;
{
  long m, k, gncell[3], stride[3], ip[2], ldim[3];
  long idim;
  FluidWork *fw;
  real *Radius, *Colatitude;
  real radius=0.0, central_radius=0.0, colatitude=0.0;
  fw = CurrentFluidPatch;
  Radius = fw->desc->Center[_RAD_];
  Colatitude = fw->desc->Center[_COLAT_];
  ldim[0] = dim;
  ldim[1] = ip[0] = (dim == 0);
  ldim[2] = ip[1] = 2-(dim == 2);
  getgridsize (fw->desc, gncell, stride);
  for (k = 1; k < gncell[dim]; k++) {
    m = ii*stride[ip[0]]+jj*stride[ip[1]]+k*stride[dim];
      fw->Flux_mass[dim][m] = beam->mass_flux[k];
    if (!Isothermal)
      fw->Flux_energy[dim][m] = beam->energy_flux[k];
    fw->InterfacePressure[dim][m] = beam->pressure_godunov[k];
    for (idim = 0; idim < NDIM; idim++){
      fw->Flux_mom[ldim[idim]][dim][m] = beam->momentum_flux[idim][k];
//if(fw->desc->level==7&&m==59354&&CPU_Rank==3&&ldim[idim]==1) printf("BEAM flux mom=%lg oth side=%lg dim=%ld ldim(idm)=%ld idim=%ld\n",beam->momentum_flux[0][k],beam->momentum_flux[0][k+stride[dim]],dim,ldim[idim],idim);
    }
				/* The following lines substitute the linear
				   momentum flux for the angular
				   momentum flux if appropriate (i.e.,
				   in cylindrical or spherical
				   coordinates, for azimuthal
				   momentum */
    if (__CYLINDRICAL) {

      radius = Radius[m];
      if (dim == _RAD_)
      	radius = fw->desc->Edges[_RAD_][k];
      if ((_AZIM_ < NDIM) && (dim != _AZIM_))
	fw->Flux_mom[_AZIM_][dim][m] = \
	  fw->Flux_mom[_AZIM_][dim][m] * radius +\
	  (beam->mass_flux[k]) * OMEGAFRAME * radius *radius;
  }
    if (__SPHERICAL) {
      central_radius = Radius[m];
      if (dim == _RAD_)
      	central_radius = fw->desc->Edges[_RAD_][k];
      colatitude = Colatitude[m];
      if (dim == _COLAT_)
      	colatitude = fw->desc->Edges[_COLAT_][k];
      radius = central_radius * sin(colatitude);
      if ((_AZIM_ < NDIM) && (dim != _AZIM_))
	fw->Flux_mom[_AZIM_][dim][m] = \
	  fw->Flux_mom[_AZIM_][dim][m] * radius +\
	  (beam->mass_flux[k]) * OMEGAFRAME * radius *radius;
      if (_COLAT_ < NDIM)
	fw->Flux_mom[_COLAT_][dim][m] *= central_radius;
  }
  }
}



void FillDiffFluxes  (dim, ii, jj, beam)
long dim, ii, jj;
Beam *beam;
{
long m, k, gncell[3], stride[3], ip[2], ldim[3];
long idim;
FluidWork *fw;
real *Radius, *Colatitude;
real radius=0.0, central_radius=0.0, colatitude=0.0;
fw = CurrentFluidPatch;
Radius = fw->desc->Center[_RAD_];
Colatitude = fw->desc->Center[_COLAT_];
ldim[0] = dim;
ldim[1] = ip[0] = (dim == 0);
ldim[2] = ip[1] = 2-(dim == 2);
getgridsize (fw->desc, gncell, stride);
for (k = 1; k < gncell[dim]; k++) {
m = ii*stride[ip[0]]+jj*stride[ip[1]]+k*stride[dim];
 fw->Flux_diff[dim][m] = beam->diff_flux[k];


}
}

void FillDiffFluxes2  (dim, ii, jj, beam)
long dim, ii, jj;
Beam *beam;
{
long m, k, gncell[3], stride[3], ip[2], ldim[3];
long idim;
FluidWork *fw;
real *Radius, *Colatitude;
real radius=0.0, central_radius=0.0, colatitude=0.0;
fw = CurrentFluidPatch;
Radius = fw->desc->Center[_RAD_];
Colatitude = fw->desc->Center[_COLAT_];
ldim[0] = dim;
ldim[1] = ip[0] = (dim == 0);
ldim[2] = ip[1] = 2-(dim == 2);
getgridsize (fw->desc, gncell, stride);
for (k = 1; k < gncell[dim]; k++) {
m = ii*stride[ip[0]]+jj*stride[ip[1]]+k*stride[dim];
 fw->Flux_diff[dim][m] = beam->diff_flux[k];


for (idim = 0; idim < NDIM; idim++){
 fw->Flux_diff_mom[ldim[idim]][dim][m] = beam->momentum_diff_flux[idim][k];

}
   /* The following lines substitute the linear
      momentum flux for the angular
      momentum flux if appropriate (i.e.,
      in cylindrical or spherical
      coordinates, for azimuthal
      momentum */
if (__CYLINDRICAL) {

 radius = Radius[m];
 if (dim == _RAD_)
   radius = fw->desc->Edges[_RAD_][k];
 if ((_AZIM_ < NDIM) && (dim != _AZIM_))
fw->Flux_diff_mom[_AZIM_][dim][m] = \
fw->Flux_diff_mom[_AZIM_][dim][m] * radius +\
(beam->diff_flux[k]) * OMEGAFRAME * radius *radius;
}
if (__SPHERICAL) {
 central_radius = Radius[m];
 if (dim == _RAD_)
   central_radius = fw->desc->Edges[_RAD_][k];
 colatitude = Colatitude[m];
 if (dim == _COLAT_)
   colatitude = fw->desc->Edges[_COLAT_][k];
 radius = central_radius * sin(colatitude);
 if ((_AZIM_ < NDIM) && (dim != _AZIM_))
fw->Flux_diff_mom[_AZIM_][dim][m] = \
fw->Flux_diff_mom[_AZIM_][dim][m] * radius +\
(beam->diff_flux[k]) * OMEGAFRAME * radius *radius;
 if (_COLAT_ < NDIM)
fw->Flux_diff_mom[_COLAT_][dim][m] *= central_radius;
}
}
}
