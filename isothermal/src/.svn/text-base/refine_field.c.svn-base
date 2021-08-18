#include "jupiter.h"

void refine_field (radix, nboutput, nbgrid, grids, gtr, cpugpatch, ncomp, cpatch, cpar)
     char radix[MAXLINELENGTH];
GridFileInfo grids[MAXGRIDS];
long nboutput, nbgrid, gtr, cpugpatch, ncomp;
tGrid_CPU *cpatch, *cpar;
{
  char srcfile[MAXLINELENGTH];
  char destfile[MAXLINELENGTH];
  FILE *in, *out;
  long chunk_length, i, j, k, l, ip, jp, kp, lp, c, nc, resolp[3], resol[3], offset[3];
  real chunk_dest[MAXLENGTHONEDIM], chunk_src[MAXLENGTHONEDIM];
  long vartype;
  long ratio, totalsize, totalsizep, size[3], sizep[3], sizeg[3], sizegp[3], lpng, lng;
  real radius_src, radius_dest, colat_src, colat_dest, interm;
  if ((_AZIM_)  && KEPLERIAN) {
    prs_stderr ("Keplerian disk on the spot refinement implemented only");
    prs_error  ("if the azimuth is the first dimension. Use COORDPERMUT accordingly");
    /* This is meant to have the same 'radius' and 'colatitude' in the
       chunks considered below */
  }
  vartype = _other_;
  if (__CYLINDRICAL)
    colat_src = colat_dest = M_PI/2.0;
  setout  (nboutput);
  sprintf (srcfile,  "%s%ld_%ld_%ld.dat", radix, nboutput, nbgrid, grids[gtr].level);
  sprintf (destfile, "%s%ld_%ld_%ld.dat", radix, nboutput, cpugpatch, grids[0].level);
  in  = prs_openrd (srcfile);
  out = prs_opend  (destfile);
  for (i = 0; i < 3; i++) {
    size[i]   = grids[0].size[i];
    sizeg[i]  = size[i]+2*Nghost[i];
    sizep[i]  = grids[gtr].size[i];
    sizegp[i] = sizep[i]+2*Nghost[i];
    resolp[i] = (grids[gtr].nc_max[i]-grids[gtr].nc_min[i])/sizep[i];
    resol[i]  = (grids[0].nc_max[i]-grids[0].nc_min[i])/size[i];
    offset[i] = (grids[0].nc_min[i]-grids[gtr].nc_min[i])/resolp[i];
  }
  ratio = (Refine[0] ? 2 : 1);
  chunk_length = size[0]*resol[0]/resolp[0];
  totalsize = size[0]*size[1]*size[2];
  totalsizep= sizep[0]*sizep[1]*sizep[2];
  for (nc = 0; nc < ncomp; nc++) {
    if (nc > 0) vartype = _other_; /* this trick can be used since if KEPLERIAN is set */
				/* then the azimuth is the first dimension */
    for (lng = 0; lng < totalsize; lng+=size[0]) {
      k = lng/(size[0]*size[1]);
      j = (lng-k*size[0]*size[1])/size[0];
      i = lng-k*size[0]*size[1]-j*size[0];
      l = (i+Nghost[0])+(j+Nghost[1])*sizeg[0]+(k+Nghost[2])*sizeg[0]*sizeg[1];
      radius_dest = cpatch->Center[_RAD_][l];
      colat_dest = cpatch->Center[_COLAT_][l];
      kp = k/(Refine[2] ? 2:1)+offset[2];
      jp = j/(Refine[1] ? 2:1)+offset[1];
      ip = i/(Refine[0] ? 2:1)+offset[0];
      lpng = ip+jp*sizep[0]+kp*sizep[0]*sizep[1];
      lp = (ip+Nghost[0])+(jp+Nghost[1])*sizegp[0]+(kp+Nghost[2])*sizegp[0]*sizegp[1];
      radius_src = cpar->Center[_RAD_][lp];
      colat_src = cpar->Center[_COLAT_][lp];
      lpng +=  totalsizep*nc;
      fseek (in, lpng*sizeof(real), SEEK_SET);
      fread (chunk_src, sizeof(real), (size_t)(chunk_length), in);
      for (c = 0; c < chunk_length*ratio; c++)
	if (KEPLERIAN) {
	  interm = keplerian_comm(chunk_src[c/ratio],vartype,radius_src,colat_src,1);
	  chunk_dest[c] = keplerian_comm(interm,vartype,radius_dest,colat_dest,-1);
	} else {
	  chunk_dest[c] = chunk_src[c/ratio];
	}
      fwrite (chunk_dest, sizeof(real), (size_t)(chunk_length*ratio), out);
    }
  }
  fclose (in);
  fclose (out);
}
