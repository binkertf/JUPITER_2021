#include "jupiter.h"

void ConservativeUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim;
  long ip1, ip2, id, istop;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  //  real delta_energy, kinetic, vlin[3];
  real *inter[3], radius, *ipres[3];
  fw = CurrentFluidPatch;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);
  for (i = 0; i < NDIM; i++) {
    inter[i] = fw->desc->InterSurface[i];
    ipres[i] = fw->InterfacePressure[i];
    momentum[i] = fw->Momentum[i];
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
  }
  _radius = fw->desc->Center[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  energy = fw->Energy;
  energy_tot = fw->Energy_tot;
  Density = fw->Density;	/* In this loop the momentum[3]  array
				   is   filled; if appropriate one  of
				   these values represents the angular
				   momentum */
  /* In the non isothermal case we also fill the total energy (internal+kinetic) for each zone */
 for (idm = 0; idm < NDIM; idm++) {
    Velf = Velocity[idm];
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  //vlin[idm] = Velocity[idm][m];
	  //momentum[idm][m] = Density[m]*vlin[idm];
	  momentum[idm][m] = Density[m]*Velf[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius; // this is the angular momentum
	    //vlin[idm] = radius*(vlin[idm]+OMEGAFRAME);
	  }
	  if (__SPHERICAL && (idm == _AZIM_)) {
	    momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    momentum[idm][m]*= radius*radius;
	    //vlin[idm] = radius*(vlin[idm]+OMEGAFRAME);
	  }
	  if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
	    radius = _radius[m];
	    momentum[idm][m]*= radius*radius;
	    //vlin[idm] *= radius;
	  }
	}
      }
    }
  }


	real visc;

	visc = VISCOSITY;
	if (fw->Fluid->next != NULL){ //dust fluid
		visc=0.0;
	}



 if (VISCOSITY > 1e-16) //for testing!
    ApplyViscousStress (dt);


//	if (!Isothermal) {	/* The kinetic energy below is evaluated in an inertial frame */
//	  kinetic = 0.0;
//	  for (idm = 0; idm < NDIM; idm++)
//	    kinetic += vlin[idm]*vlin[idm];
//	  energy_tot[m] = energy[m]+Density[m]*.5*kinetic;	/* Volumic total energy */
//	}
  for (mdim = 0; mdim < NDIM; mdim++) {
      /* In this loop the momenta are updated from the fluxes given by
	 the Riemann solver, as well as the total energy in the non
	 isothermal version */
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
	  // if (mdim == _COLAT_) delta_mom /= _radius[m];
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
	// delta_energy = 0.0;
	delta_mass=0.0;
	for (idm = 0; idm < NDIM; idm++) {
	  m_other_side = m+stride[idm];
	  delta_mass+= fw->Flux_mass[idm][m];
	  delta_mass-= fw->Flux_mass[idm][m_other_side];
	//	    if (!Isothermal) {
	//	      delta_energy+= fw->Flux_energy[idm][m];
	//	      delta_energy-= fw->Flux_energy[idm][m_other_side];
	//	    }
	}
	Density[m] += delta_mass*InvVolume[m];
	//	  if (!Isothermal)
	//	    energy_tot[m] += delta_energy*InvVolume[m];
      }
    }
  }
  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
    ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	       idx[1] = j;
	        for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	           idx[0] = i;
	            m = i*stride[0]+j*stride[1]+k*stride[2];
	  //	  vlin[idm] = Velocity[idm][m] = momentum[idm][m]/Density[m];
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	    //  vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
	  } else if (__SPHERICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    Velocity[idm][m] = momentum[idm][m]/(radius*radius);
	    Velocity[idm][m] /= Density[m];
	    Velocity[idm][m] -= OMEGAFRAME;
	    //	    vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
	  } else if (__SPHERICAL && (idm == _COLAT_)) { /* Colatitude velocity is angular */
	    Velocity[idm][m] = momentum[idm][m]/Density[m];
	    radius = _radius[m];
	    Velocity[idm][m] /= (radius * radius);
	    //	    vlin[idm] = Velocity[idm][m]*radius;
	  } else
	    Velocity[idm][m] = momentum[idm][m]/Density[m];

//	if (!Isothermal) {	/* The kinetic energy below is evaluated in an inertial frame */
//	  kinetic = 0.0;
//	  for (idm = 0; idm < NDIM; idm++)
//	    kinetic += vlin[idm]*vlin[idm];
//	  energy[m] = energy_tot[m]-.5*kinetic*Density[m];
//	}
	  if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][INF][id] += fw->Flux_mass[idm][m];
	    //	    if (!Isothermal)
	    //	      fw->Fluid->EnergyFlux->Flux[idm][INF][id] += fw->Flux_energy[idm][m];
	    if (idm != _COLAT_)
	      fw->Fluid->Pressure->Pressure[idm][INF][id] += ipres[idm][m]*inter[idm][m];
	    else
	      fw->Fluid->Pressure->Pressure[idm][INF][id] += ipres[idm][m]*inter[idm][m]*_radius[m];
	    for (mdim = 0; mdim < NDIM; mdim++) {
	//if ((idm == _COLAT_) && (mdim == _COLAT_))
	//	fw->Fluid->MomentumFlux[mdim]->Flux[idm][INF][id] += fw->Flux_mom[mdim][idm][m]*_radius[m];
	//else
		fw->Fluid->MomentumFlux[mdim]->Flux[idm][INF][id] += fw->Flux_mom[mdim][idm][m];
	    }
	  }
	  if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
	    ip1 = (idm == 0);
	    ip2 = 2-(idm == 2);
	    mp = m+stride[idm];
	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
	    fw->Fluid->MassFlux->Flux[idm][SUP][id] += fw->Flux_mass[idm][mp];
	    if (!Isothermal)
	      fw->Fluid->EnergyFlux->Flux[idm][SUP][id] += fw->Flux_energy[idm][mp];
	    if (idm != _COLAT_)
	      fw->Fluid->Pressure->Pressure[idm][SUP][id] += ipres[idm][mp]*inter[idm][mp];
	    else
	      fw->Fluid->Pressure->Pressure[idm][SUP][id] += ipres[idm][mp]*inter[idm][mp]*_radius[m];
	    /* _radius[m] or mp does not make a difference in the above sentence */
	    for (mdim = 0; mdim < NDIM; mdim++) {
	      //if ((idm == _COLAT_) && (mdim == _COLAT_))
	      //fw->Fluid->MomentumFlux[mdim]->Flux[idm][INF][id] += fw->Flux_mom[mdim][idm][mp]*_radius[m];
	      /* _radius[m] or mp does not make a difference in the above sentence */
	      //else
		fw->Fluid->MomentumFlux[mdim]->Flux[idm][SUP][id] += fw->Flux_mom[mdim][idm][mp];
	    }
	  }
	}
      }
    }
  }
}






void DiffusionUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim;
  long ip1, ip2, id, istop;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf, Velo;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  real *inter[3], radius, *ipres[3];

  fw = CurrentFluidPatch;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);

  for (i = 0; i < NDIM; i++) {
    inter[i] = fw->desc->InterSurface[i];
    ipres[i] = fw->InterfacePressure[i];
    momentum[i] = fw->Momentum[i];
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
  }
  _radius = fw->desc->Center[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  Density = fw->Density; /* In this loop the densities are updated from the fluxes given by
     the diffusion solver */





  //momentum
    for (idm = 0; idm < NDIM; idm++) {
        Velf = Velocity[idm];
        for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
          for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
              m = i*stride[0]+j*stride[1]+k*stride[2];
                momentum[idm][m] = Density[m]*Velf[m];
                if (__CYLINDRICAL && (idm == _AZIM_)) {
                  momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
                  radius = _radius[m];
                  momentum[idm][m]*= radius*radius; // this is the angular momentum
                }
                if (__SPHERICAL && (idm == _AZIM_)) {
                  momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
                  radius = _radius[m];
                  radius *= sin(_colatitude[m]);
                  momentum[idm][m]*= radius*radius;
                }
                if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
                  radius = _radius[m];
                  momentum[idm][m]*= radius*radius;
                }
              }
            }
          }
        }

        for (mdim = 0; mdim < NDIM; mdim++) { //loop over the momenta
            /* In this loop the momenta are updated from the diffusion fluxes */
          for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
            for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
                for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
                   delta_mom = 0.0;
                   m = i*stride[0]+j*stride[1]+k*stride[2];
                   for (idm = 0; idm < NDIM; idm++) { //loop over the directions
                      m_other_side = m+stride[idm];
                      delta_mom+= fw->Flux_diff_mom[mdim][idm][m];
                      delta_mom-= fw->Flux_diff_mom[mdim][idm][m_other_side];
                    }
                    momentum[mdim][m] += delta_mom*InvVolume[m];
                 }
              }
            }
          }



          // mass
            for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
              for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
                for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {

                  m = i*stride[0]+j*stride[1]+k*stride[2];

                  delta_mass=0.0;
                  for (idm = 0; idm < NDIM; idm++) {
                    m_other_side = m+stride[idm];

                    //if (fw->Flux_diff[idm][m]>0.0)
                      //printf("dtg2: %lg \n",fw->Flux_diff[idm][m]);

                    delta_mass+= fw->Flux_diff[idm][m];
                    delta_mass-= fw->Flux_diff[idm][m_other_side];
                    }



                  Density[m] += delta_mass*InvVolume[m];
                  //if (Density[m]<DUSTDENSFLOOR)
                      //Density[m] = DUSTDENSFLOOR;

                }
              }
            }


  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
    ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
         idx[1] = j;
          for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
             idx[0] = i;
              m = i*stride[0]+j*stride[1]+k*stride[2];


              if (__CYLINDRICAL && (idm == _AZIM_)) {
                   radius = _radius[m];
                   Velocity[idm][m] = momentum[idm][m]/(radius*radius);
                   Velocity[idm][m] /= Density[m];
                   Velocity[idm][m] -= OMEGAFRAME;
                   //  vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
                 } else if (__SPHERICAL && (idm == _AZIM_)) {
                   radius = _radius[m];
                   radius *= sin(_colatitude[m]);
                   Velocity[idm][m] = momentum[idm][m]/(radius*radius);
                   Velocity[idm][m] /= Density[m];
                   Velocity[idm][m] -= OMEGAFRAME;
                   //      vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
                 } else if (__SPHERICAL && (idm == _COLAT_)) { /* Colatitude velocity is angular */
                   Velocity[idm][m] = momentum[idm][m]/Density[m];
                   radius = _radius[m];
                   Velocity[idm][m] /= (radius * radius);
                   //      vlin[idm] = Velocity[idm][m]*radius;
                 } else
                   Velocity[idm][m] = momentum[idm][m]/Density[m];
            }
          }
        }
      }







  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
    ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	       idx[1] = j;
	        for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	           idx[0] = i;
	            m = i*stride[0]+j*stride[1]+k*stride[2];



 if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
    fw->Fluid->DiffFlux->Flux[idm][INF][id] += fw->Flux_diff[idm][m];
    }

    if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
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


              for (idm = 0; idm < NDIM; idm++) {
                ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
                ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
                for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
                  idx[2] = k;
                  for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            	       idx[1] = j;
            	        for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
            	           idx[0] = i;
            	            m = i*stride[0]+j*stride[1]+k*stride[2];



             if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
                id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];

                for (mdim = 0; mdim < NDIM; mdim++) {
                  //fw->Fluid->MomentumDiffFlux[mdim]->Flux[idm][INF][id] += fw->Flux_diff_mom[mdim][idm][m];
                  fw->Fluid->MomentumDiffFlux[mdim]->Flux[idm][INF][id] += fw->Flux_diff_mom[idm][mdim][m];
                }
              }


              if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
            	    ip1 = (idm == 0);
            	    ip2 = 2-(idm == 2);
            	    mp = m+stride[idm];
            	    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
                  for (mdim = 0; mdim < NDIM; mdim++) {
                    //fw->Fluid->MomentumDiffFlux[mdim]->Flux[idm][SUP][id] += fw->Flux_diff_mom[mdim][idm][mp];
                    fw->Fluid->MomentumDiffFlux[mdim]->Flux[idm][SUP][id] += fw->Flux_diff_mom[idm][mdim][mp];
                  }
            	 }



            }
            }
            }
            }


}

void MassUpdate (real dt)
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, idx[3], size[3];
  long m_other_side, m, mp, idm, mdim;
  long ip1, ip2, id, istop;
  real *momentum[3], *Velocity[3], *Density, *InvVolume, *Velf, Velo;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real delta_mass, delta_mom;
  real *inter[3], radius, *ipres[3];

  fw = CurrentFluidPatch;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);

  for (i = 0; i < NDIM; i++) {
    inter[i] = fw->desc->InterSurface[i];
    ipres[i] = fw->InterfacePressure[i];
    momentum[i] = fw->Momentum[i];
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
  }
  _radius = fw->desc->Center[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  Density = fw->Density; /* In this loop the densities are updated from the fluxes given by
     the diffusion solver */





  //momentum
    for (idm = 0; idm < NDIM; idm++) {
        Velf = Velocity[idm];
        for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
          for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
              m = i*stride[0]+j*stride[1]+k*stride[2];
                momentum[idm][m] = Density[m]*Velf[m];
                if (__CYLINDRICAL && (idm == _AZIM_)) {
                  momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
                  radius = _radius[m];
                  momentum[idm][m]*= radius*radius; // this is the angular momentum
                }
                if (__SPHERICAL && (idm == _AZIM_)) {
                  momentum[idm][m] = Density[m]*(Velf[m]+OMEGAFRAME);
                  radius = _radius[m];
                  radius *= sin(_colatitude[m]);
                  momentum[idm][m]*= radius*radius;
                }
                if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
                  radius = _radius[m];
                  momentum[idm][m]*= radius*radius;
                }
              }
            }
          }
        }


  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
    ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
         idx[1] = j;
          for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
             idx[0] = i;
              m = i*stride[0]+j*stride[1]+k*stride[2];


              if (__CYLINDRICAL && (idm == _AZIM_)) {
                   radius = _radius[m];
                   Velocity[idm][m] = momentum[idm][m]/(radius*radius);
                   Velocity[idm][m] /= Density[m];
                   Velocity[idm][m] -= OMEGAFRAME;
                   //  vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
                 } else if (__SPHERICAL && (idm == _AZIM_)) {
                   radius = _radius[m];
                   radius *= sin(_colatitude[m]);
                   Velocity[idm][m] = momentum[idm][m]/(radius*radius);
                   Velocity[idm][m] /= Density[m];
                   Velocity[idm][m] -= OMEGAFRAME;
                   //      vlin[idm] = (Velocity[idm][m]+OMEGAFRAME)*radius;
                 } else if (__SPHERICAL && (idm == _COLAT_)) { /* Colatitude velocity is angular */
                   Velocity[idm][m] = momentum[idm][m]/Density[m];
                   radius = _radius[m];
                   Velocity[idm][m] /= (radius * radius);
                   //      vlin[idm] = Velocity[idm][m]*radius;
                 } else
                   Velocity[idm][m] = momentum[idm][m]/Density[m];
            }
          }
        }
      }


      // mass
        for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
          for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {

              m = i*stride[0]+j*stride[1]+k*stride[2];

              delta_mass=0.0;
              for (idm = 0; idm < NDIM; idm++) {
                m_other_side = m+stride[idm];

                //if (fw->Flux_diff[idm][m]>0.0)
                  //printf("dtg2: %lg \n",fw->Flux_diff[idm][m]);

                delta_mass+= fw->Flux_diff[idm][m];
                delta_mass-= fw->Flux_diff[idm][m_other_side];
                }



              Density[m] += delta_mass*InvVolume[m];
              if (Density[m]<DUSTDENSFLOOR)
                  Density[m] = DUSTDENSFLOOR;

            }
          }
        }




  for (idm = 0; idm < NDIM; idm++) {
    ip1 = (idm == 0);	/* If yes then keep track of flux and interface pressure */
    ip2 = 2-(idm == 2);	/* in order to ensure conservation across nested grids */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      idx[2] = k;
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	       idx[1] = j;
	        for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	           idx[0] = i;
	            m = i*stride[0]+j*stride[1]+k*stride[2];



 if (idx[idm] == Nghost[idm]) {	/* Are we on the lower face ? */
    id  = (idx[ip1]-Nghost[ip1])+(idx[ip2]-Nghost[ip2])*size[ip1];
    fw->Fluid->DiffFlux->Flux[idm][INF][id] += fw->Flux_diff[idm][m];
    }

    if (idx[idm] == size[idm]+Nghost[idm]-1) { /* Are we on the upper face ? */
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
