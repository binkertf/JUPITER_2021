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
    __Riemann_Solver (rhoL, rhoR, uL, uR, aL, aR, &us, &ps);
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
    if (( uL < 0.0 ) && ( uR > 0.0 )) { /* the two states are moving away from each other leaving vacuum in between. The flux is zero */
      for (k = 0; k < NDIM-1; k++){
        vp_i[k] = 0.0;
        }
      rho_i = 0.0; //1.0e-20;
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
    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)
      beam->momentum_flux[k+1][i] = rho_i*vp_i[k]*u_i * surfdt;
    beam->pressure_godunov[i] = 0.0;
  }
}


real sgn(real x){
  if (x >= 0.0) return 1.0;
  if (x < 0.0) return -1.0;
}

real wavelimit(real s[3], real Z[3])
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
    beam->mass_flux[i] = mass_flux * surfdt;
    beam->momentum_flux[0][i] = pu_flux * surfdt;
    for (k = 0; k < NDIM-1; k++)
      beam->momentum_flux[k+1][i] = pup_flux[k] * surfdt;
    beam->pressure_godunov[i] = 0.0;
  }
}


void Compute_Fluxes_Diffusion(beam,beam2,dt)
     Beam *beam, *beam2;
     real dt;
{
  long i,dim,k;
  real rho_i, pu_flux, mass_flux, uf, radius, am;
  real rhoL, rhoR, rhogL, rhogR, dtgL, dtgR,uL, uR, u_i, vp_i[2], vpL[2], vpR[2];
  real D_d,Pi,surfdt, j, dx, omegakep, St, cs, v_d;
  dim = beam->dim[0];
  for (i = Nghost[dim]; i <= beam->length-Nghost[dim]; i++) {
    rhoL = beam->rhoL[i]; //dust density
    rhoR = beam->rhoR[i]; //dust density

    rhogL = beam2->rhoL[i]; //gas density
    rhogR = beam2->rhoR[i]; //gas density

    dx = beam->center[i]-beam->center[i-1];
    D_d = VISCOSITY;


    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _RAD_)) { //account for epicycle oscillations in radial direction
      if (constSt!=TRUE){
        radius = beam->center[i];
        omegakep = 1.0/(sqrt(radius*radius*radius)); //local keplerian frequeny (at midplane)OMEGAFRAME/(sqrt(radius*radius*radius)); //local keplerian frequeny (at midplane)
        St = sqrt(M_PI/8.0)*DUSTSIZE*omegakep/(sqrt(rhogL*rhogR)*radius*ASPECTRATIO); //average Stokes number
      }else{
        St = STOKESNUMBER;
      }
      D_d = D_d/(1.0+St*St); // see Youdin&Lithwick (2007)
    }

    //if ((__CYLINDRICAL || __SPHERICAL) && (dim == _COLAT_)) { //account for epicycle oscillations in radial direction
      //St = STOKESNUMBER;
      //D_d = D_d/(1.0+St*St); // see Youdin&Lithwick (2007)
    //}

    Pi = -D_d*sqrt((rhoL+rhogL)*(rhoR+rhogR))*(rhoR/(rhogR+rhoR)-rhoL/(rhogL+rhoL))/dx; //diffusion flux (geometric mean)



    // Diffusion flux limiter
    cs = sqrt(beam2->cs[i-1]*beam2->cs[i]); //gas sound speed
    v_d = sqrt(Pi*Pi)/sqrt(rhoL*rhoR); //diffusion velocity
    if(v_d>cs){
      Pi = Pi/sqrt(Pi*Pi) * cs*sqrt(rhoL*rhoR);
    }



    if ((__CYLINDRICAL || __SPHERICAL) && (dim == _AZIM_)) { //add geometrical correction
    Pi = Pi / beam->radius;
    }

    surfdt = beam->intersurface[i] * dt;
    beam->diff_flux[i] = Pi * surfdt; //dust mass diffusion flux

  }

}
