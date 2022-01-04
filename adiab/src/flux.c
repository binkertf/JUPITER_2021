#include "jupiter.h"

void ConservativeUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim, len;
  long mrp, mcp;
  long ip1, ip2, id, lev;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  real delta_energy, vlin;
  real delta_energy_tot;
  real *inter[3], radius, *ipres[3];
  real *rhov[3], *centrif, *coriolis, *radint;
  fw = CurrentFluidPatch;
  lev = fw->desc->level;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);
  len = 1;
  for (i = 0; i < NDIM; i++) {
    inter[i] = fw->desc->InterSurface[i];
    ipres[i] = fw->InterfacePressure[i];
    momentum[i] = fw->Momentum[i];
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
    len *= gncell[i];
  }
  _radius = fw->desc->Center[_RAD_];
  radint  = fw->desc->Edges[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  energy = fw->Energy;
  energy_tot = fw->Energy_tot;
  Density = fw->Density;	/* In this loop the momentum[3]  array
				   is   filled; if appropriate one  of
				   these values represents the angular
				   momentum */
  /* In the non isothermal case we also fill the total energy (internal+kinetic) for each zone */
  centrif  = fw->Centrifugal_acc;
  coriolis = fw->Coriolis;
  for (idm = 0; idm < NDIM; idm++) {
    Velf = Velocity[idm];
    rhov[idm] = fw->RawMassFlux[idm];
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  vlin = Velocity[idm][m]; /* Linear velocity *in corotating frame* */
	  momentum[idm][m] = Density[m]*Velf[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius;
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    momentum[idm][m]*= radius*radius;
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius;
	    vlin *= radius;
	  }
	}
      }
    }
  }

  /* The function below must be AFTER the calculation of momenta,
     because the geometric terms alter directly the value of the
     momenta. */
  if (VISCOSITY > 1e-16)
    JUP_SAFE(ApplyViscousStress (dt));
  /* The function below must be AFTER the flux update due to the
     viscosity, as the viscous stress alters the momenta fluxes. It
     must also be placed, obviously, BEFORE the conservative update
     performed in the following loop. */
  JUP_SAFE(CorrectFluxesFromFinerLevel ());
  /* the fluxes are then corrected from finer level */

  for (mdim = 0; mdim < NDIM; mdim++) {
    /* In this loop the momenta are updated from the fluxes given by
       the Riemann solver*/
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  delta_mom = 0.0;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  for (idm = 0; idm < NDIM; idm++) {
	    m_other_side = m+stride[idm];
	    delta_mom+= fw->Flux_mom[mdim][idm][m];
	    delta_mom-= fw->Flux_mom[mdim][idm][m_other_side];
	  }
	  momentum[mdim][m] += delta_mom*InvVolume[m];
	}
      }
    }
  }
  /* In this loop the densities are updated from the fluxes given by
     the  Riemann solver, as  well as  the total  energy in  the non
     isothermal version */
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	delta_energy_tot = 0.0;
	delta_energy     = 0.0;
	delta_mass=0.0;
	for (idm = 0; idm < NDIM; idm++) {
	  m_other_side = m+stride[idm];
	  delta_mass+= fw->Flux_mass[idm][m];
	  delta_mass-= fw->Flux_mass[idm][m_other_side];
	  if (!Isothermal) {
	    delta_energy+= fw->Flux_energy[idm][m];
	    delta_energy-= fw->Flux_energy[idm][m_other_side];
	    delta_energy_tot += fw->Flux_tot_energy[idm][m];
	    delta_energy_tot -= fw->Flux_tot_energy[idm][m_other_side];
	  }
	}
	Density[m] += delta_mass*InvVolume[m];
	if (!Isothermal) {
	  mrp = m + stride[_RAD_];
	  if (_COLAT_ < NDIM)
	    mcp = m + stride[_COLAT_];
	  energy_tot[m] += delta_energy_tot*InvVolume[m];
	  if (lev < HIGHRESLEVEL) {
	    if (__CYLINDRICAL) {
	      energy_tot[m] -= .5*(rhov[_RAD_][m]  *centrif[m]+		\
				   rhov[_RAD_][mrp]*centrif[mrp])*dt;
	      energy_tot[m] -= (coriolis[m]+coriolis[mrp])*dt*OMEGAFRAME;
	    }
	    if (__SPHERICAL) {
	      energy_tot[m] -= .5*(rhov[_RAD_][m]*centrif[m]+		\
				   rhov[_RAD_][mrp]*centrif[mrp])*dt;
	      energy_tot[m] -= (coriolis[m]+coriolis[mrp])*dt*OMEGAFRAME;
	      energy_tot[m] += .5*(rhov[_COLAT_][m]+rhov[_COLAT_][mcp])*OMEGAFRAME*OMEGAFRAME* \
		(_radius[m])*sin(_colatitude[m])*cos(_colatitude[m])*dt;
	    }
	  }
	}
      }
    }
  }
  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);
    ip2 = 2-(idm == 2);
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	idx[1] = j;
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  idx[0] = i;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  Velocity[idm][m] = momentum[idm][m]/Density[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	  } else if (__SPHERICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	  } else if (__SPHERICAL && (idm == _COLAT_)) { /* Colatitude velocity is angular */
	    Velocity[idm][m] = momentum[idm][m]/Density[m];
	    radius = _radius[m];
	    Velocity[idm][m] /= (radius * radius);
	  }
	  if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
	    /* If yes then keep track of flux and interface pressure
	       in order to ensure conservation across nested grids */
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][INF][id] += fw->Flux_mass[idm][m];
	    if (!Isothermal) {
	      fw->Fluid->TotalEnergyFlux->Flux[idm][INF][id] += fw->Flux_tot_energy[idm][m];
	    }
	    fw->Fluid->Pressure->Pressure[idm][INF][id] += ipres[idm][m]*inter[idm][m];
	    for (mdim = 0; mdim < NDIM; mdim++) {
	      fw->Fluid->MomentumFlux[mdim]->Flux[idm][INF][id] += fw->Flux_mom[mdim][idm][m];
	    }
	  }
	  if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
	    /* If yes then keep track of flux and interface pressure
	       in order to ensure conservation across nested grids */
	    ip1 = (idm == 0);
	    ip2 = 2-(idm == 2);
	    mp = m+stride[idm];
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][SUP][id] += fw->Flux_mass[idm][mp];
	    if (!Isothermal) {
	      fw->Fluid->TotalEnergyFlux->Flux[idm][SUP][id] += fw->Flux_tot_energy[idm][mp];
	    }
	    fw->Fluid->Pressure->Pressure[idm][SUP][id] += ipres[idm][mp]*inter[idm][mp];
	    for (mdim = 0; mdim < NDIM; mdim++) {
	      fw->Fluid->MomentumFlux[mdim]->Flux[idm][SUP][id] += fw->Flux_mom[mdim][idm][mp];
	    }
	  }
	}
      }
    }
  }
}




void ConservativeDustUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim, len;
  long mrp, mcp;
  long ip1, ip2, id, lev;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  real delta_energy, vlin;
  real delta_energy_tot;
  real *inter[3], radius, *ipres[3];
  real *rhov[3], *centrif, *coriolis, *radint;
  fw = CurrentFluidPatch;
  lev = fw->desc->level;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);
  len = 1;
  for (i = 0; i < NDIM; i++) {
    inter[i] = fw->desc->InterSurface[i];
    ipres[i] = fw->InterfacePressure[i];
    momentum[i] = fw->Momentum[i];
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
    len *= gncell[i];
  }
  _radius = fw->desc->Center[_RAD_];
  radint  = fw->desc->Edges[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  Density = fw->Density;	/* In this loop the momentum[3]  array
				   is   filled; if appropriate one  of
				   these values represents the angular
				   momentum */

  for (idm = 0; idm < NDIM; idm++) {
    Velf = Velocity[idm];
    rhov[idm] = fw->RawMassFlux[idm];
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  vlin = Velocity[idm][m]; /* Linear velocity *in corotating frame* */
	  momentum[idm][m] = Density[m]*Velf[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius;
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    momentum[idm][m]*= radius*radius;
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius;
	    vlin *= radius;
	  }
	}
      }
    }
  }

  /* The function below must be AFTER the calculation of momenta,
     because the geometric terms alter directly the value of the
     momenta. */
  //if (VISCOSITY > 1e-16)
    //JUP_SAFE(ApplyViscousStress (dt));

  /* The function below must be AFTER the flux update due to the
     viscosity, as the viscous stress alters the momenta fluxes. It
     must also be placed, obviously, BEFORE the conservative update
     performed in the following loop. */
  JUP_SAFE(CorrectFluxesFromFinerLevel ());
  /* the fluxes are then corrected from finer level */

  for (mdim = 0; mdim < NDIM; mdim++) {
    /* In this loop the momenta are updated from the fluxes given by
       the Riemann solver*/
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  delta_mom = 0.0;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  for (idm = 0; idm < NDIM; idm++) {
	    m_other_side = m+stride[idm];
	    delta_mom+= fw->Flux_mom[mdim][idm][m];
	    delta_mom-= fw->Flux_mom[mdim][idm][m_other_side];
	  }
	  momentum[mdim][m] += delta_mom*InvVolume[m];
	}
      }
    }
  }
  /* In this loop the densities are updated from the fluxes given by
     the  Riemann solver, as  well as  the total  energy in  the non
     isothermal version */
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	delta_mass=0.0;
	for (idm = 0; idm < NDIM; idm++) {
	  m_other_side = m+stride[idm];
	  delta_mass+= fw->Flux_mass[idm][m];
	  delta_mass-= fw->Flux_mass[idm][m_other_side];

	}
	Density[m] += delta_mass*InvVolume[m];

      }
    }
  }
  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);
    ip2 = 2-(idm == 2);
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	idx[1] = j;
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  idx[0] = i;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  Velocity[idm][m] = momentum[idm][m]/Density[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	  } else if (__SPHERICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	  } else if (__SPHERICAL && (idm == _COLAT_)) { /* Colatitude velocity is angular */
	    Velocity[idm][m] = momentum[idm][m]/Density[m];
	    radius = _radius[m];
	    Velocity[idm][m] /= (radius * radius);
	  }
	  if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
	    /* If yes then keep track of flux and interface pressure
	       in order to ensure conservation across nested grids */
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][INF][id] += fw->Flux_mass[idm][m];
	    //if (!Isothermal) {
	      //fw->Fluid->TotalEnergyFlux->Flux[idm][INF][id] += fw->Flux_tot_energy[idm][m];
	    //}
	    //fw->Fluid->Pressure->Pressure[idm][INF][id] += ipres[idm][m]*inter[idm][m];
	    for (mdim = 0; mdim < NDIM; mdim++) {
	      fw->Fluid->MomentumFlux[mdim]->Flux[idm][INF][id] += fw->Flux_mom[mdim][idm][m];
	    }
	  }
	  if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
	    /* If yes then keep track of flux and interface pressure
	       in order to ensure conservation across nested grids */
	    ip1 = (idm == 0);
	    ip2 = 2-(idm == 2);
	    mp = m+stride[idm];
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][SUP][id] += fw->Flux_mass[idm][mp];
	    //if (!Isothermal) {
	      //fw->Fluid->TotalEnergyFlux->Flux[idm][SUP][id] += fw->Flux_tot_energy[idm][mp];
	    //}
	    //fw->Fluid->Pressure->Pressure[idm][SUP][id] += ipres[idm][mp]*inter[idm][mp];
	    for (mdim = 0; mdim < NDIM; mdim++) {
	      fw->Fluid->MomentumFlux[mdim]->Flux[idm][SUP][id] += fw->Flux_mom[mdim][idm][mp];
	    }
	  }
	}
      }
    }
  }
}





void  DiffusionUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim, len;
  long mrp, mcp;
  long ip1, ip2, id, lev;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  real delta_energy, vlin;
  real delta_energy_tot;
  real *inter[3], radius, *ipres[3];
  real *rhov[3], *centrif, *coriolis, *radint;
  fw = CurrentFluidPatch;
  lev = fw->desc->level;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);
  len = 1;
  for (i = 0; i < NDIM; i++) {
    size[i] = gncell[i]-2*Nghost[i];
    //size[i] = gncell[i];
    len *= gncell[i];
  }
  _radius = fw->desc->Center[_RAD_];
  radint  = fw->desc->Edges[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];

  Density = fw->Density;



  /* In this loop the densities are updated from the fluxes given by
     the diffusion solver */
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	         m = i*stride[0]+j*stride[1]+k*stride[2];

	         delta_mass=0.0;
	          for (idm = 0; idm < NDIM; idm++) {
	             m_other_side = m+stride[idm];
	             delta_mass+= fw->Flux_diff[idm][m];
	             delta_mass-= fw->Flux_diff[idm][m_other_side];
	           }
	           Density[m] += delta_mass*InvVolume[m];
             if (Density[m]<DUSTDENSFLOOR){
               Density[m] = DUSTDENSFLOOR;
             }

      }
    }
  }

  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);
    ip2 = 2-(idm == 2);
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	       idx[1] = j;
         for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
            idx[0] = i;
            m = i*stride[0]+j*stride[1]+k*stride[2];

	          if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
	                /* If yes then keep track of flux and interface pressure
	                 in order to ensure conservation across nested grids */
	                  id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
                    fw->Fluid->DiffFlux->Flux[idm][INF][id] += fw->Flux_diff[idm][m];
            }
            if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
	          /* If yes then keep track of flux and interface pressure
	          in order to ensure conservation across nested grids */
              ip1 = (idm == 0);
	            ip2 = 2-(idm == 2);
	            mp = m+stride[idm];
	            id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
              fw->Fluid->DiffFlux->Flux[idm][SUP][id] += fw->Flux_diff[idm][mp];

	     }
	    }
      }
      }
      }
}
