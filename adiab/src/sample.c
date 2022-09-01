#include "jupiter.h"

void Compute_Fluxes_Iso (beam, dt)
     Beam *beam;
     real dt;
{
  long i,dim,k;
  real rhoL, rhoR, aL, aR, uL, uR, vpL[2], vpR[2];
  real us, ps, S_shock, a_godunov;
  real rho_i, u_i, vp_i[2]; 	/* _i like 'interface' */
  real mass_flux, pu_flux;
  real surfdt;
  real a2, uf, v0=0.0, radius;
  /* Azimuthal flux corrections for Keplerian disks below */
  real C_mass_1=1.0, C_mass_2=1.0;
  real C_mom_1=1.0, C_mom_2=1.0;
  boolean HSE;
  dim = beam->dim[0];
  HSE = CorrHydroStat[dim];
  if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
    radius = beam->radius;
    v0 = radius * OMEGAFRAME;
    if (KEPLERIAN && MERIDIANCORRECTION) {
      C_mass_1 = beam->masscorr1;
      C_mass_2 = beam->masscorr2;
      C_mom_1 = beam->momcorr1;
      C_mom_2 = beam->momcorr2;
    }
  }
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    uL   = beam->uL[i];
    uR   = beam->uR[i];
    for (k = 0; k < NDIM-1; k++) {
      vpL[k] = beam->v_perp_L[k][i];
      vpR[k] = beam->v_perp_R[k][i];
    }
    rhoL = beam->rhoL[i];
    rhoR = beam->rhoR[i];
    if (HSE) {
      a2 = beam->cs2i[i];
      aR = aL   = sqrt(a2);
    } else {
      aR = aL   = .5*(beam->cs[i-1]+beam->cs[i]);
      a2 = aL * aL;
    }
    /* Call the Riemann solver */
    GetStar_TWOSHOCKS(rhoL, rhoR, uL, uR, aL, aR, &us, &ps);
    if (0.0 < us) {
      /* we are in the star or left region or fan */
      for (k = 0; k < NDIM-1; k++)
	vp_i[k] = vpL[k];
      a_godunov = aL;
      if (us < uL) {		/* There is a left shock */
	S_shock = uL-sqrt(ps/rhoL);
	if (0.0 < S_shock) {	/* We are at the left of the left shock */
	  u_i = uL;
	  rho_i = rhoL;
	} else {		/* we are in the star region */
	  u_i = us;
	  rho_i = ps/a2;
	}
      } else {			/* There is a left rarefaction */
	if (0.0 < uL-aL) {	/* We are in the left region */
	  u_i = uL;
	  rho_i = rhoL;
	} else {
	  if (0.0 > us-aL) {	/* We are in the star region */
	    u_i = us;
	    rho_i = ps/a2;
	  } else {		/* We are in the fan */
	    u_i = aL;
	    rho_i = rhoL*exp((uL-u_i)/aL);
	  }
	}
      }
    } else {			/* We are in the star or right region or fan */
      for (k = 0; k < NDIM-1; k++)
	vp_i[k] = vpR[k];
      a_godunov = aR;
      if (us > uR) {		/* There is a right shock */
	S_shock = uR+sqrt(ps/rhoR);
	if (0.0 > S_shock) {	/* we are at the right of the right shock */
	  u_i = uR;
	  rho_i = rhoR;
	} else {		/* we are in the star region */
	  u_i = us;
	  rho_i = ps/a2;
	}
      } else {			/* There is a right rarefaction */
	if (0.0 > uR+aR) {	/* we are in the right part */
	  u_i = uR;
	  rho_i = rhoR;
	} else {
	  if (0.0 < us+aR) {	/* we are in the star region */
	    u_i = us;
	    rho_i = ps/a2;
	  } else {		/* we are in the fan */
	    u_i = -aR;
	    rho_i = rhoR*exp(-(uR-u_i)/aR);
	  }
	}
      }
    }
    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
      uf = u_i+v0;
      mass_flux =rho_i*(uf*C_mass_1-v0*C_mass_2);
      pu_flux = (a2+uf*uf)*C_mom_1-v0*uf*C_mom_2;
      pu_flux *= radius * rho_i;
    } else {
      mass_flux =rho_i*u_i;
      pu_flux   = rho_i*(a2+u_i*u_i);
    }
    surfdt = beam->intersurface[i] * dt;
    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)
      beam->momentum_flux[k+1][i] = rho_i*vp_i[k]*u_i * surfdt;
    beam->pressure_godunov[i] = rho_i*a2;
  }
}

void Compute_Fluxes_Adi (beam, dt)
     Beam *beam;
     real dt;
{
  long i,dim,ip[2],k;
  real rhoL, rhoR, eL, eR, aL, aR, uL, uR, pL, pR, vpL[2], vpR[2], gammaL, gammaR;
  real us, ps, S_shock, rhoS, aS, gammaS;
  real rho_i, u_i, p_i, a_i, vp_i[2], gamma_i; 	// _i like 'interface'
  real ekp; // perpendicular kinetic energy (over rho)
  real mass_flux, pu_flux, ene_flux;
  real surfdt, SLStar, SRStar;
  real uf, v0=0.0, radius, wratio;
  // Azimuthal flux corrections for Keplerian disks below
  real C_mass_1=1.0, C_mass_2=1.0;
  real C_mom_1=1.0, C_mom_2=1.0;
  real C_ene_1=1.0, C_ene_2=1.0;
  boolean VacuumCreated;
  dim = beam->dim[0];
  ip[0] = beam->dim[1];
  ip[1] = beam->dim[2];
  if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
    radius = beam->radius;
    v0 = radius * OMEGAFRAME; // Multiplication by sin(colatitude) already included in beam.c, for the spherical case.
    if (KEPLERIAN && MERIDIANCORRECTION) {
      C_mass_1 = beam->masscorr1;
      C_mass_2 = beam->masscorr2;
      C_mom_1 = beam->momcorr1;
      C_mom_2 = beam->momcorr2;
      C_ene_1 = beam->enercorr1;
      C_ene_2 = beam->enercorr2;
    }
  }
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    uL   = beam->uL[i];
    uR   = beam->uR[i];
    for (k = 0; k < NDIM-1; k++) {
      vpL[k] = beam->v_perp_L[k][i];
      vpR[k] = beam->v_perp_R[k][i];
    }
    rhoL = beam->rhoL[i];
    rhoR = beam->rhoR[i];
    eL = beam->eL[i];
    eR = beam->eR[i];

    if(!Isothermal && Stellar && !CONST_GAMMA){
      pL = secant_method(eL, rhoL, 0.00001); // pressure
      pR = secant_method(eR, rhoR, 0.00001); // pressure
      gammaL = Gamma1(pL * MU / rhoL * TEMP0, rhoL);
      gammaR = Gamma1(pR * MU / rhoR * TEMP0, rhoR);
    }
    
    else{
      gammaL = GAMMA;
      gammaR = GAMMA;
      pL = beam->eL[i]*(gammaL-1.); // pressure
      pR = beam->eR[i]*(gammaR-1.); // pressure
    }    
    
    aR = sqrt(gammaL*pR/rhoR); // sound speed
    aL = sqrt(gammaR*pL/rhoL); // sound speed
    // Call the Riemann solver
    VacuumCreated =   GetStar_AdiabaticSolver (rhoL, rhoR, uL, uR, aL, aR, &us, &ps, gammaL, gammaR);
    if (VacuumCreated) {
      SRStar = uR - aR * 2.0/(gammaR-1.0); 
      SLStar = uL + aL * 2.0/(gammaL-1.0);
      // Since vacuum is created, we know that SRStar > SLStar. Vacuum
      // resides in between these two characteristics. On each side,
      // we have rarefaction waves.  See Toro, p.142, fig. 4.17 and
      // the corresponding paragraph.  We sample the solution fan in
      // this case.
      if (SLStar > 0.0) { //We are in the left region or left fan:
			  //vacuum is in right space
	      for (k = 0; k < NDIM-1; k++)
	        vp_i[k] = vpL[k];
	      if (uL - aL > 0.0) { // Left state
	        u_i = uL;
	        rho_i = rhoL;
	        a_i = aL;
          gamma_i = gammaL;
	      }
        else {  // Left fan
	        wratio = (gammaL-1.0)/(gammaL+1.0)*uL/aL;
	        u_i = 2.0/(gammaL+1.0) * (aL + uL*(gammaL-1.)*0.5);
	        rho_i = rhoL * pow(2.0/(gammaL+1.0)  + wratio, 2.0/(gammaL-1.0));
	        p_i = pL * pow(2.0/(gammaL+1.0) + wratio, gammaL*2.0/(gammaL-1.0));
	        a_i = sqrt(gammaL*p_i/rho_i);
          gamma_i = gammaL;
	      }
      }
      if (SRStar < 0.0) { //We are in the right region or right fan:
	                  //vacuum is in left space.  This condition
	                  //is mutually exclusive of the previous one.
	      for (k = 0; k < NDIM-1; k++)
	        vp_i[k] = vpR[k];
	      if (uR + aR < 0.0) { // Right state
	        u_i = uR;
	        rho_i = rhoR;
	        a_i = aR;
          gamma_i = gammaR;
	      }
        else { // Right fan
	        wratio = (gammaR-1.0)/(gammaR+1.0)*uR/aR;
	        u_i = 2.0/(gammaR+1.0) * (-aR + uR*(gammaR-1.)*0.5);
	        rho_i = rhoR * pow(2.0/(gammaR+1.0) - wratio, 2.0/(gammaR-1.0));
	        p_i = pR * pow(2.0/(gammaR+1.0) - wratio, gammaR*2.0/(gammaR-1.0));
	        a_i = sqrt(gammaR*p_i/rho_i);
          gamma_i = gammaR;
	      }
      }
      if ((SLStar <= 0.0) && (SRStar >= 0.0)) {// Vacuum is on the interface.
	      for (k = 0; k < NDIM-1; k++)
	        vp_i[k] = .5*(vpL[k]+vpR[k]); // This value does not matter: associated fluxes vanish
	      u_i = 0.0; // This is the linear interpolation of SLStar and
		    // SRStar at the interface :)
	      rho_i = 0.0;
	      p_i = a_i = 0.0;
      }
    }
    else {
      // Schematic code : I=Interface ; S=Shock ; D=contact Discontinuity ; F=Fan
      if (0.0 < us) { //0.0, as we need the solution on the interface
	      // we are in the star or left region or fan
	      for (k = 0; k < NDIM-1; k++)
	        vp_i[k] = vpL[k];
	      if (ps > pL) {            // There is a left shock
	        S_shock = uL + (pL - ps)/((uL - us)*rhoL);
	        wratio = ps/pL;
	        rhoS = rhoL * (wratio + (gammaL-1.0)/(gammaL+1.0))  // rhoS(on shock side)
	        /(wratio * (gammaL-1.0)/(gammaL+1.0) + 1.);
	        aS = sqrt(gammaL*ps/rhoS);
          gammaS = gammaL;
	        if (0.0 < S_shock) {    // We are at the left of the left shock: ISD
	          u_i = uL;
	          rho_i = rhoL;
	          a_i = aL;
            gamma_i = gammaL;
	        }
          else {                // we are in the left star region: SID
	          u_i = us;
	          rho_i = rhoS;
	          a_i = aS;
            gamma_i = gammaS;
	        }
	      }
        else {                    // There is a left rarefaction
	        rhoS = rhoL * pow(ps/pL, 1.0/gammaL); // rhoS(on rarefaction side)
	        aS = sqrt(gammaL*ps/rhoS);
          gammaS = gammaL;
	        if (0.0 < uL-aL) {      // We are in the left region: IFD
	          u_i = uL;
	          rho_i = rhoL;
	          a_i = aL;
            gamma_i = gammaL;
	        }
          else {
	          if (0.0 > us-aS) {  // We are in the left star region FID
	            u_i = us;
	            rho_i = rhoS;
	            a_i = aS;
              gamma_i = gammaS;
	          }
            else {            // We are in the fan: FIFD
	            wratio = (gammaL-1.0)/(gammaL+1.0)*uL/aL;
	            u_i = 2.0/(gammaL+1.0) * (aL + uL*(gammaL-1.)*0.5);
	            rho_i = rhoL * pow(2.0/(gammaL+1.0) + wratio, 2.0/(gammaL-1.0));
	            p_i = pL * pow(2.0/(gammaL+1.0) + wratio, gammaL*2.0/(gammaL-1.0));
	            a_i = sqrt(gammaL*p_i/rho_i);
              gamma_i = gammaL;
	          }
	        }
	      }
      } 
      else {                        // We are in the star or right region or fan
	      for (k = 0; k < NDIM-1; k++)
	        vp_i[k] = vpR[k];
	      if (ps > pR) {            // There is a right shock
	        S_shock = uR + (pR - ps)/(uR - us)/rhoR;
	        wratio = ps/pR;
	        rhoS = rhoR * (wratio + (gammaR-1.0)/(gammaR+1.0)) //rhoS (on shock side)
	        /(wratio * (gammaR-1.0)/(gammaR+1.0) + 1.);
	        aS = sqrt(gammaR*ps/rhoS);
          gammaS = gammaR;
	        if (0.0 > S_shock) {    // we are at the right of the right shock: DSI
	          u_i = uR;
	          rho_i = rhoR;
	          a_i = aR;
            gamma_i = gammaR;
	        }
          else {                // we are in the right star region: DIS
	          u_i = us;
	          rho_i = rhoS;
	          a_i = aS;
            gamma_i = gammaS;
	        }
	      }
        else {                    // There is a right rarefaction
	        rhoS = rhoR * pow(ps/pR,1.0/gammaR); //rhoS(on rarefaction side)
	        aS = sqrt(gammaR*ps/rhoS);
          gammaS = gammaR;
	        if (0.0 > uR+aR) {      // we are in the right part: DFI
	          u_i = uR;
	          rho_i = rhoR;
	          a_i = aR;
            gamma_i = gammaR;
	        }
          else {
	          if (0.0 < us+aS) {  // we are in the right star region: DIF
	            u_i = us;
	            rho_i = rhoS;
	            a_i = aS;
              gamma_i = gammaS;
	          }
            else {            // we are in the fan: DFIF
	            wratio = (gammaR-1.0)/(gammaR+1.0)*uR/aR;
	            u_i = 2.0/(gammaR+1.0) * (-aR + uR*(gammaR-1.)*0.5);
	            rho_i = rhoR * pow(2.0/(gammaR+1.0) - wratio, 2.0/(gammaR-1.0));
	            p_i = pR * pow(2.0/(gammaR+1.0) - wratio, gammaR*2.0/(gammaR-1.0));
	            a_i = sqrt(gammaR*p_i/rho_i);
              gamma_i = gammaR;
	          }
	        }
	      }
      }
    }
    // Modified by Judit: we set the colat & radial velocties to zero in the first active cells:
    if((dim == _COLAT_ && i==2 && beam->rawcoord[i] < RANGE3LOW+(beam->rawcoord[1]-beam->rawcoord[0])) || (dim == _RAD_ && i==2 && beam->rawcoord[i] < RANGE2LOW+.5*(beam->rawcoord[1]-beam->rawcoord[0]))){
      u_i = 0.0; // CHECK IF OK (beambound.c)
    }
    // modification made till here
    ekp = 0.5*u_i*u_i;
    for (k = 0; k < NDIM-1; k++)
      ekp += vp_i[k]*vp_i[k]*0.5;
    if (dim == _RAD_) {
      beam->centrifugal_acc[i] = vp_i[0]*vp_i[0];
      if (_COLAT_ < NDIM)
	      beam->centrifugal_acc[i] += vp_i[1]*vp_i[1];
      beam->centrifugal_acc[i] /= beam->edge[i];
      beam->coriolis[i] = rho_i*vp_i[0]*u_i;
      if (__SPHERICAL)
	      beam->coriolis[i] *= sin(beam->colatitude);
    }
    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
      uf = u_i+v0;
      mass_flux = rho_i * (uf*C_mass_1 - v0*C_mass_2);
      pu_flux   = (a_i*a_i/gamma_i + uf*uf)*C_mom_1 - v0*uf*C_mom_2;
      pu_flux  *= radius * rho_i;     // angular momentum in inertial frame...
      ene_flux = uf*rho_i * (a_i*a_i/(gamma_i-1.0) + ekp)* C_ene_1-	\
	    v0*rho_i * (a_i*a_i/(gamma_i-1.0) + ekp)*C_ene_2;
    }
    else {
      mass_flux = rho_i * u_i;
      pu_flux   = rho_i * (a_i*a_i/gamma_i + u_i*u_i);
      ene_flux  = u_i*rho_i * (a_i*a_i/(gamma_i-1.0) + ekp);
    }

    surfdt = beam->intersurface[i] * dt;
    beam->raw_mass_flux[i] = mass_flux;
    beam->u_interface[0][i] = u_i;
    if (NDIM > 1)
      beam->u_interface[1][i] = vp_i[0];
    if (NDIM > 2)
      beam->u_interface[2][i] = vp_i[1];
    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)    // will be corrected in FillFluxes
      beam->momentum_flux[k+1][i] = rho_i*vp_i[k]*u_i * surfdt;
    beam->tot_energy_flux[i] = ene_flux * surfdt;
    beam->pressure_godunov[i] = a_i*a_i*rho_i*gamma_i;
  }
}



void Compute_Fluxes_pressureless1 (beam, dt) /* fisrt order pressureless flux */
   Beam *beam;
     real dt;
{
  long i,dim,k;
  real rhoL, rhoR, aL, aR, uL, uR, vpL[2], vpR[2];
  real us, ps, S_shock, a_godunov;
  real rho_i, u_i, vp_i[2]; 	/* _i like 'interface' */
  real mass_flux, pu_flux;
  real surfdt;
  real a2, uf, v0=0.0, radius;
  /* Azimuthal flux corrections for Keplerian disks below */
  real C_mass_1=1.0, C_mass_2=1.0;
  real C_mom_1=1.0, C_mom_2=1.0;
  boolean HSE;
  dim = beam->dim[0];
  HSE = CorrHydroStat[dim];
  if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
    radius = beam->radius;
    v0 = radius * OMEGAFRAME; /* linear keplerian velcoity v0 */
    if (KEPLERIAN && MERIDIANCORRECTION) {
      C_mass_1 = beam->masscorr1;
      C_mass_2 = beam->masscorr2;
      C_mom_1 = beam->momcorr1;
      C_mom_2 = beam->momcorr2;
    }
  }
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    uL   = beam->uL[i];
    uR   = beam->uR[i];
    for (k = 0; k < NDIM-1; k++) {
      vpL[k] = beam->v_perp_L[k][i];
      vpR[k] = beam->v_perp_R[k][i];
    }
    rhoL = beam->rhoL[i];
    rhoR = beam->rhoR[i];

    aR = aL = a2 = 0.0; /* sound speed is zero for pressurless fluid */
    /* Call the Riemann solver */
    GetStar_PRESSURELESS(rhoL, rhoR, uL, uR, aL, aR, &us, &ps); /* the pressureless solver calculates the delta wave velocity uL */
    if (( uL < 0.0 ) && ( uR > 0.0 )){ /* the two states are moving away from each other leaving vacuum in between. The flux is zero */
      for (k = 0; k < NDIM-1; k++){
        vp_i[k] = 0.0;
        }
      rho_i = 0.0;//rho_i = 1.0e-20;
      u_i = 0.0;

    } else {
      if ( us > 0.0) { /* Delta wave is moving to the right*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = vpL[k];
        }
        rho_i = rhoL;
        u_i = uL;
      }
      if ( us < 0.0) { /* Delta wave is moving to the left*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = vpR[k];
        }
        rho_i = rhoR;
        u_i = uR;
      }
      if ( us == 0.0){ /* Delta wave stays at the interface*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = 0.5*(vpL[k]+vpR[k]);
        }
        rho_i = 0.5*(rhoR+rhoL);
        u_i = 0.5*(uL+uR);
      }
    }
 
if(!Isothermal){
          if (dim == _RAD_) {
      beam->centrifugal_acc[i] = vp_i[0]*vp_i[0];
      if (_COLAT_ < NDIM)
	beam->centrifugal_acc[i] += vp_i[1]*vp_i[1];
      beam->centrifugal_acc[i] /= beam->edge[i];
      beam->coriolis[i] = rho_i*vp_i[0]*u_i;
      if (__SPHERICAL)
	beam->coriolis[i] *= sin(beam->colatitude);
    }
}

      if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) { /* compute the flux vector G here (azimuthal flux) and save as angular flux */
      uf = u_i+v0; /* linear azimuthal velocity*/
      mass_flux =rho_i*(uf*C_mass_1-v0*C_mass_2); /*density flux */
      pu_flux = (uf*uf)*C_mom_1-v0*uf*C_mom_2; /*momentum flux (pressureless)*/
      pu_flux *= radius * rho_i;



    } else { /* compute the flux vector F here (radial flux) */
      mass_flux = rho_i*u_i; /* density flux*/
      pu_flux   = rho_i*(u_i*u_i); /* momentum flux (pressureless)*/
    }
    surfdt = beam->intersurface[i] * dt;

    beam->raw_mass_flux[i] = mass_flux;
    beam->u_interface[0][i] = u_i;
    if (NDIM > 1)
      beam->u_interface[1][i] = vp_i[0];
    if (NDIM > 2)
      beam->u_interface[2][i] = vp_i[1];

    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)
      beam->momentum_flux[k+1][i] = rho_i*vp_i[k]*u_i * surfdt;
      if (!Isothermal){
        beam->tot_energy_flux[i] = 0.0;
      }
      beam->pressure_godunov[i] = 0.0;
  }
}


// ################################################################

real sgn(real x){
  if (x >= 0.0) return 1.0;
  if (x < 0.0) return -1.0;
}

real wavelimit(real s[3], real Z[3]) //minmod
{
   real Zl,Zr,Zc,r,a,phi;
    Zc = 0.0;
    Zl = Z[0];
    Zr = Z[2];
    Zc = Z[1];

    if (s[1]>0.0){
      r = Zl/Zc;
    }else{
      r = Zr/Zc;
    }

    phi = 0.0;

    if (r<1.0){
      phi = r;
    } else{
      phi = 1.0;
    }
    if (r < 0.0){
      phi = 0.0;
    }

  return Zc * phi;
}


void Compute_Fluxes_pressureless2 (beam, dt) /* higher order pressureless flux */
     Beam *beam;
     real dt;
{
  long i,dim,k,r;
  real rhoL, rhoR, aL, aR, uL, uR, vpL[2], vpR[2];
  real us, ps, S_shock, a_godunov;
  real rho_i, u_i, vp_i[2]; 	/* _i like 'interface' */
  real mass_flux, pu_flux, pup_flux[2];
  real cmass_flux, cpu_flux, cpup_flux[2]; /* correction fluxes */
  real fL_mass, fL_pu, fR_mass, fR_pu, fL_pup[2], fR_pup[2];
  real Z1_mass[3], Z1_pu[3], Z1_pup[3][2],Z2_mass[3], Z2_pu[3], Z2_pup[3][2],xpup[3];
  real s1[3],s2[3],dx;
  real surfdt;
  real a2, uf, v0=0.0, radius;
  /* Azimuthal flux corrections for Keplerian disks below */
  real C_mass_1=1.0, C_mass_2=1.0;
  real C_mom_1=1.0, C_mom_2=1.0;
  boolean HSE;
  dim = beam->dim[0];
  HSE = CorrHydroStat[dim];
  if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
    radius = beam->radius;
    v0 = radius * OMEGAFRAME;
    if (KEPLERIAN && MERIDIANCORRECTION) {
      C_mass_1 = beam->masscorr1;
      C_mass_2 = beam->masscorr2;
      C_mom_1 = beam->momcorr1;
      C_mom_2 = beam->momcorr2;
    }
  }
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    uL   = beam->uL[i];
    uR   = beam->uR[i];
    for (k = 0; k < NDIM-1; k++) {
      vpL[k] = beam->v_perp_L[k][i];
      vpR[k] = beam->v_perp_R[k][i];
    }
    rhoL = beam->rhoL[i];
    rhoR = beam->rhoR[i];

    aR = aL = a2 = 0.0; /* sound speed is zero for pressurless fluid */

    /* Call the Riemann solver */
    GetStar_PRESSURELESS(rhoL, rhoR, uL, uR, aL, aR, &us, &ps); /* the pressureless solver calculates the delta wave velocity uL */
    if (( uL < 0.0 ) && ( uR > 0.0 )) { /* the two states are moving away from each other leaving vacuum in between. The flux is zero */
      for (k = 0; k < NDIM-1; k++){
        vp_i[k] = 0.0;
        }
      rho_i = 0.0;
      u_i = 0.0;

    } else {
      if ( us > 0.0) { /* Delta wave is moving to the right*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = vpL[k];
        }
        rho_i = rhoL;
        u_i = uL;
      }
      if ( us < 0.0) { /* Delta wave is moving to the left*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = vpR[k];
        }
        rho_i = rhoR;
        u_i = uR;
      }
      if ( us == 0.0){ /* Delta wave stays at the interface*/
        for (k = 0; k < NDIM-1; k++){
          vp_i[k] = 0.5*(vpL[k]+vpR[k]);
        }
        rho_i = 0.5*(rhoR+rhoL);
        u_i = 0.5*(uL+uR);
      }
    }

    /*----------------------------------------------------------------*/
    /* compute correction fluxes */
    /*----------------------------------------------------------------*/

  for (r=0;r<=2;r++){

    uL   = beam->uL[i-1+r];
    uR   = beam->uR[i-1+r];
    for (k = 0; k < NDIM-1; k++) {
      vpL[k] = beam->v_perp_L[k][i-1+r];
      vpR[k] = beam->v_perp_R[k][i-1+r];
    }
    rhoL = beam->rhoL[i-1+r];
    rhoR = beam->rhoR[i-1+r];

    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
      fL_mass =rhoL*((uL+v0)*C_mass_1-v0*C_mass_2);
      fL_pu = rhoL*radius*((uL+v0)*(uL+v0)*C_mom_1-v0*(uL+v0)*C_mom_2);
      fR_mass = rhoR*((uR+v0)*C_mass_1-v0*C_mass_2);
      fR_pu = rhoR*radius*((uR+v0)*(uR+v0)*C_mom_1-v0*(uR+v0)*C_mom_2);
    } else {
      fL_mass =rhoL*uL;
      fL_pu = rhoL*uL*uL;
      fR_mass = rhoR*uR;
      fR_pu = rhoR*uR*uR;
    }


    for (k = 0; k < NDIM-1; k++){
          fL_pup[k]=rhoL*uL*vpL[k];
          fR_pup[k]=rhoR*uR*vpR[k];
        }

    /* Call the Riemann solver */
    GetStar_PRESSURELESS(rhoL, rhoR, uL, uR, aL, aR, &us, &ps); /* the pressureless solver calculates the delta wave velocity uL */
    if (( uL < 0.0 ) && ( uR > 0.0 )) { /* the two states are moving away from each other leaving vacuum in between. The flux is zero */
      for (k = 0; k < NDIM-1; k++){
        Z1_pup[r][k]=-fL_pup[k];
        Z2_pup[r][k]=fR_pup[k];
        }
      Z1_mass[r] = -fL_mass;
      Z1_pu[r] = -fL_pu;
      Z2_mass[r] = fR_mass;
      Z2_pu[r] = fR_pu;
      s1[r] = uL;
      s2[r] = uR;
    } else {
      if ( us > 0.0) { /* Delta wave is moving to the right*/
        for (k = 0; k < NDIM-1; k++){
          Z2_pup[r][k]=fR_pup[k]-fL_pup[k];
          Z1_pup[r][k]=0.0;
        }
        Z2_mass[r] = fR_mass-fL_mass;
        Z2_pu[r] = fR_pu-fL_pu;
        Z1_mass[r] = 0.0;
        Z1_pu[r] = 0.0;
        s1[r] = us;
        s2[r] = us;

      }
      if ( us < 0.0) { /* Delta wave is moving to the left*/
        for (k = 0; k < NDIM-1; k++){
        Z1_pup[r][k]=fR_pup[k]-fL_pup[k];
        Z2_pup[r][k]=0.0;
        }
        Z1_mass[r] = fR_mass-fL_mass;
        Z1_pu[r] = fR_pu-fL_pu;
        Z2_mass[r] = 0.0;
        Z2_pu[r] = 0.0;
        s1[r] = us;
        s2[r] = us;
      }
      if ( us == 0.0){ /* Delta wave stays at the interface*/
        for (k = 0; k < NDIM-1; k++){
          Z2_pup[r][k]=fR_pup[k]-fL_pup[k];
          Z1_pup[r][k]=0.0;
        }
        Z2_mass[r] = fR_mass-fL_mass;
        Z2_pu[r] = fR_pu-fL_pu;
        Z1_mass[r] = 0.0;
        Z1_pu[r] = 0.0;
        s1[r] = us;
        s2[r] = us;
      }
    }

    }

      Z1_mass[1] = wavelimit(s1,Z1_mass);
      Z2_mass[1] = wavelimit(s2,Z2_mass);
      Z1_pu[1] = wavelimit(s1,Z1_pu);
      Z2_pu[1] = wavelimit(s2,Z2_pu);
      for (k = 0; k < NDIM-1; k++){
          xpup[0]= Z1_pup[0][k];
          xpup[1]= Z1_pup[1][k];
          xpup[2]= Z1_pup[2][k];
          Z1_pup[1][k] = wavelimit(s1,xpup);
          xpup[0]= Z2_pup[0][k];
          xpup[1]= Z2_pup[1][k];
          xpup[2]= Z2_pup[2][k];
          Z2_pup[1][k] = wavelimit(s2,xpup);
        }


    dx = beam->edge[i+1]-beam->edge[i];
    for (k = 0; k < NDIM-1; k++){
          cpup_flux[k] = (0.5*sgn(s1[1])*(1.0-dt/dx*sqrt(s1[1]*s1[1]))*Z1_pup[1][k])+(0.5*sgn(s2[1])*(1.0-dt/dx*sqrt(s2[1]*s2[1]))*Z2_pup[1][k]);
        }
    cmass_flux = (0.5*sgn(s1[1])*(1.0-dt/dx*sqrt(s1[1]*s1[1]))*Z1_mass[1])+(0.5*sgn(s2[1])*(1.0-dt/dx*sqrt(s2[1]*s2[1]))*Z2_mass[1]);
    cpu_flux = (0.5*sgn(s1[1])*(1.0-dt/dx*sqrt(s1[1]*s1[1]))*Z1_pu[1])+(0.5*sgn(s2[1])*(1.0-dt/dx*sqrt(s2[1]*s2[1]))*Z2_pu[1]);



/*----------------------------------------------------------------------*/


if(!Isothermal){

          if (dim == _RAD_) {
      beam->centrifugal_acc[i] = vp_i[0]*vp_i[0];
      if (_COLAT_ < NDIM)
	beam->centrifugal_acc[i] += vp_i[1]*vp_i[1];
      beam->centrifugal_acc[i] /= beam->edge[i];
      beam->coriolis[i] = rho_i*vp_i[0]*u_i;
      if (__SPHERICAL)
	beam->coriolis[i] *= sin(beam->colatitude);
    }
}
    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) {
      uf = u_i+v0;
      mass_flux =rho_i*(uf*C_mass_1-v0*C_mass_2); /*density flux*/
      pu_flux = (uf*uf)*C_mom_1-v0*uf*C_mom_2; /*momentum flux (pressureless)*/
      pu_flux *= radius * rho_i;


      mass_flux+=cmass_flux;
      pu_flux+=(cpu_flux);

    } else {
      mass_flux =rho_i*u_i; /* density flux*/
      pu_flux   = rho_i*(u_i*u_i); /* momentum flux (pressureless)*/

      mass_flux = mass_flux + cmass_flux;
      pu_flux = pu_flux + cpu_flux;
    }




    for (k = 0; k < NDIM-1; k++){
          pup_flux[k] = (rho_i*vp_i[k]*u_i) + cpup_flux[k];
        }

    surfdt = beam->intersurface[i] * dt;

    beam->raw_mass_flux[i] = mass_flux;
    beam->u_interface[0][i] = u_i;
    if (NDIM > 1)
      beam->u_interface[1][i] = vp_i[0];
    if (NDIM > 2)
      beam->u_interface[2][i] = vp_i[1];


    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)
      beam->momentum_flux[k+1][i] = pup_flux[k] * surfdt;
      if (!Isothermal){
        beam->tot_energy_flux[i] = 0.0;
      }
      beam->pressure_godunov[i] = 0.0;
  }
}





void Compute_Fluxes_Diffusion(beam,beam2,dt)
     Beam *beam, *beam2;
     real dt;
{
  long i,dim,k;
  real rhoL, rhoR, rhogL, rhogR, dtgL, dtgR, gammaL, gammaR, eL, eR, pL, pR;
  real D_d,Pi,surfdt, j, dx, radius, omegakep, St, cs, csL, csR, v_d, v_d_max, t_diff, t_stop, dustsz, dustsolidrho;
  
  dustsz = DUSTSIZE / R0; // diameter of dust grains in code units
  dustsolidrho = DUSTSOLIDRHO / RHO0; // solid density of dust grains in code units, typically 3 g/cm^3 in physical units

  dim = beam->dim[0];
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    rhoL = beam->rhoL[i]; //dust density
    rhoR = beam->rhoR[i]; //dust density

    rhogL = beam2->rho[i-1]; //gas density
    rhogR = beam2->rho[i]; //gas density

    eL = beam->eL[i];
    eR = beam->eR[i];

    dx = beam->center[i]-beam->center[i-1];
    D_d = VISCOSITY;


    if(!Isothermal && Stellar && !CONST_GAMMA){
      pL = secant_method(eL, rhoL, 0.00001); // pressure
      pR = secant_method(eR, rhoR, 0.00001); // pressure
      gammaL = Gamma1(pL * MU / rhoL * TEMP0, rhoL);
      gammaR = Gamma1(pR * MU / rhoR * TEMP0, rhoR);
      printf("gammaL: %f, gammaR: %f \n", gammaL, gammaR);
    }
    
    else{
      gammaL = GAMMA;
      gammaR = GAMMA;
      pL = beam->eL[i]*(gammaL-1.); // pressure
      pR = beam->eR[i]*(gammaR-1.); // pressure
    }


    if (!Isothermal){
      csL = sqrt(beam2->cs[i-1]/rhogL*gammaL*(gammaL-1.0)); //adiabatic sound speed
      csR = sqrt(beam2->cs[i]/rhogR*gammaR*(gammaR-1.0)); //adiabatic sound speed
      if(!constSt==YES){
        t_stop = sqrt(M_PI*sqrt(gammaL*gammaR)/8.0)*dustsolidrho*dustsz/(sqrt(rhogL*rhogR)*sqrt(csL*csR));; //local stopping time
      }
    }else{
      csL = beam2->cs[i-1];
      csR = beam2->cs[i];
      if(!constSt==YES){
        t_stop = sqrt(M_PI/8.0)*dustsolidrho*dustsz/(sqrt(rhogL*rhogR)*sqrt(csL*csR));; //local stopping time
      }
    }

    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _RAD_)) { //account for epicycle oscillations in radial direction
      if(constSt==YES){
        St = STOKESNUMBER;
      }
      else{
        radius = beam->center[i];
        omegakep = 1.0/(sqrt(radius)*sqrt(radius)*sqrt(radius)); //local keplerian frequeny (at midplane)
        St = t_stop * omegakep; //local Stokes number
      }
      D_d = D_d/(1.0+St*St); // see Youdin&Lithwick (2007)
      
    }

    //diffusion time-scale limit 1 - this is likely too conservative
    //t_diff = dx * dx / D_d;  //diffusioon timescale
    //if (t_diff < t_stop){
    //  D_d = dx * dx / t_stop;
    //}


    Pi = -D_d*sqrt((rhoL+rhogL)*(rhoR+rhogR))*(rhoR/(rhogR+rhoR)-rhoL/(rhogL+rhoL))/dx; //diffusion flux
    //Pi = -D_d*sqrt((rhogL)*(rhogR))*(rhoR/(rhogR)-rhoL/(rhogL))/dx; //diffusion flux

    // Diffusion flux limiter - the diffusion velocity must be smaller than the sound speed in gas
    cs = 0.1 * sqrt(csL * csR);
    v_d = sqrt(Pi * Pi) / sqrt(rhoL * rhoR); //diffusion velocity
    if(v_d > cs){
      Pi = sgn(Pi) * cs * sqrt(rhoL * rhoR);
    }

    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) { //add geometrical correction
      Pi = Pi / beam->radius;
    }

    surfdt = beam->intersurface[i] * dt;
    beam->diff_flux[i] = Pi * surfdt; //dust diffusion flux
  }

}
