#include "jupiter.h"

real BB_Xmin[3], BB_Xmax[3];

void Ray_System_BoundingBox () {
  long j;
  if (__CARTESIAN) {
    for (j = 0; j < 3; j++) {
      BB_Xmin[j] = corner_min0[j];
      BB_Xmax[j] = corner_max0[j];
    }
  }
  if (__CYLINDRICAL) {
    BB_Xmin[0] = -corner_max0[_RAD_];
    BB_Xmax[0] = +corner_max0[_RAD_];
    BB_Xmin[1] = -corner_max0[_RAD_];
    BB_Xmax[1] = +corner_max0[_RAD_];
    BB_Xmin[2] = corner_min0[_Vertical_];
    BB_Xmax[2] = +corner_max0[_Vertical_];
  }
  if (__SPHERICAL) {
    BB_Xmin[0] = -corner_max0[_RAD_];
    BB_Xmax[0] = +corner_max0[_RAD_];
    BB_Xmin[1] = -corner_max0[_RAD_];
    BB_Xmax[1] = +corner_max0[_RAD_];
    BB_Xmin[2] = corner_max0[_RAD_]*cos(corner_max0[_COLAT_]);
    BB_Xmax[2] = corner_max0[_RAD_]*cos(corner_min0[_COLAT_]);
    if (RayMirror)
      BB_Xmin[2] = -BB_Xmax[2];
  }
  pInfo ("System bounding box:\n");
  pInfo ("Xmin, Ymin, Zmin = %g, %g, %g\n",\
	 BB_Xmin[0], BB_Xmin[1], BB_Xmin[2]);
  pInfo ("Xmax, Ymax, Zmax = %g, %g, %g\n",\
	 BB_Xmax[0], BB_Xmax[1], BB_Xmax[2]);
}

void Ray_Tracing () {
  long output, count=0, Rx, Ry, depth, channels;
  char camera_file[MAXLINELENGTH];
  char name[MAXLINELENGTH], mirror[MAXLINELENGTH];
  char s[30][MAXLINELENGTH], st[MAXLINELENGTH];
  FILE *cam, *pic;
  real Xo, Yo, Zo, Xt, Yt, Zt, Cx, Cy, Ox, Oy, Oz;
  real ux, uy, uz, u;
  real vx, vy, vz;
  real nx, ny, nz, n;
  real dx, dy, dz, C1, C2, res;
  float *data_cube, *ray;
  long i, j ;
  if (CPU_Number > 1) {
    prs_error ("Ray tracing must be run on one processor only.\n");
  }
  if (NDIM < 3) {
    prs_error ("NDIM must be 3 for ray tracing...\n");
  } 
  sscanf (RayTracingInfo, "%ld %s", &output, camera_file);
  cam = prs_openr (camera_file, "r");
  while ((fgets(st, 199, cam) != NULL) && (count < 30)) {
    if (st[0] != '#') {
      strncpy (s[count], st, MAXLINELENGTH-1);
      count++;
    }
  }
  sscanf(s[0], "%lf %lf %lf", &Xo, &Yo, &Zo); /* o like observer */
  sscanf(s[1], "%lf %lf %lf", &Xt, &Yt, &Zt); /* t like target */
  sscanf(s[2], "%ld %ld", &Rx, &Ry); /* Camera resolution in pixels */
  sscanf(s[3], "%lf %lf", &Cx, &Cy); /* Camera aperture half angle */
  sscanf(s[4], "%lf %lf %lf", &Ox, &Oy, &Oz); /* Orientation angle */
  sscanf(s[5], "%s", name); /* Picture name */
  sscanf(s[6], "%ld", &channels); /* Number of opacity channels */
  depth = channels+9;
  sprintf (mirror, "no");
  if (count > 7)
    sscanf(s[7], "%s", mirror); /* Equatorial mirror */
  if (CoarseRayTracing) {
    Rx = 106;
    Ry = 80;
  }
  mirror[0] = (char)(toupper((int)mirror[0]));
  if ((mirror[0]=='Y') || (mirror[0]=='T') || (mirror[0]=='1'))
    RayMirror = TRUE;
  Ray_System_BoundingBox ();
  prs_msg ("Observer is located at (%g, %g, %g)\n", Xo, Yo, Zo);
  prs_msg ("Target is located at (%g, %g, %g)\n", Xt, Yt, Zt);
  prs_msg ("Camera resolution is %ld x %ld\n", Rx, Ry);
  prs_msg ("Camera aperture angle is %g x %g\n", Cx, Cy);
  prs_msg ("Camera orientation vector is (%g,%g,%g)\n", Ox, Oy, Oz);
  if (RayMirror)
    prs_msg ("Equatorial mirroring activated\n");
  else
    prs_msg ("Equatorial mirroring desactivated\n");
  prs_msg ("%ld opacity channels, %ld layers in total\n", channels, depth);
  res = sqrt((4.*Cx*Cx/(real)(Rx*Rx))+(4.*Cy*Cy/(real)(Ry*Ry)));
  data_cube = prs_malloc (Rx*Ry*sizeof(float)*depth);
  /* (ux, uy, uz) is the vector of unit length going from observer */
  /* towards the target : vec u = vec(OT)/OT */
  ux = Xt-Xo;
  uy = Yt-Yo;
  uz = Zt-Zo;
  u = sqrt(ux*ux+uy*uy+uz*uz);
  ux /= u;
  uy /= u;
  uz /= u;
  nx = uy*Oz-uz*Oy;
  ny = uz*Ox-ux*Oz;
  nz = ux*Oy-uy*Ox;
  n = sqrt(nx*nx+ny*ny+nz*nz);
  if (fabs(n) < 1e-10) {
    prs_error ("Fatal degeneracy: camera orientation vector and line of sight are aligned.\n");
  }
  nx /= n;
  ny /= n;
  nz /= n;
  vx = -uy*nz+uz*ny;
  vy = -uz*nx+ux*nz;
  vz = -ux*ny+uy*nx;
  ResetChronometer (3);
  ray = data_cube;
  for (j = 0; j < Ry; j++) {
    C2 = (-1.0+2.0*(real)j/(real)Ry)*Cy;
    for (i = 0; i < Rx; i++) {
      C1 = (-1.0+2.0*(real)i/(real)Rx)*Cx;
      dx = ux + nx*C1 + vx*C2;
      dy = uy + ny*C1 + vy*C2;
      dz = uz + nz*C1 + vz*C2;
      Ray_Integrate (Xo, Yo, Zo, dx, dy, dz, ray, depth, res);
      ray += depth;
    }
  }
  pic = prs_open (name);
  fwrite (data_cube, Rx*Ry*depth, sizeof(float), pic);
  fclose (pic);
  free (data_cube);
  ReadChronometer (3, "one picture construction");
  fclose (cam);
  prs_msg ("Written file %s, with %ld channels of size %ldx%ld\n",\
	   name, depth, Rx, Ry);
  strcat (name, ".txt");
  pic = prs_open (name);
  fprintf (pic, "%ld %ld %ld\n", depth, Rx, Ry);
  fclose (pic);
  prs_exit (1);
}
