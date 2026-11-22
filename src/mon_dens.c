//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void mon_dens_cpu () {

//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Slope);
  INPUT(Mdot);
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* interm = Slope->field_cpu;
  real* mdot = Mdot->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real gammasink = GAMMASINK;
  real rsink = RSINK;
  real dsink = DELTASINK;
  real omega_bin = BBHOMEGARATIO + OMEGAFRAME0;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real sink;
  real dx;
  real dy;
  real dz=0.0;
  real dist;
//<\INTERNAL>

//<CONSTANT>
// real Xplanet(1);
// real Yplanet(1);
// real Zplanet(1);
// real VXplanet(1);
// real VYplanet(1);
// real VZplanet(1);
// real MplanetVirtual(1);
// real Syk(Nz+2*NGHZ);
// real InvVj(Ny+2*NGHY);
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real Sxi(Nx);
//<\CONSTANT>


//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z-1; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
  ll = l;  
  interm[ll] = 0.;

#ifdef CYLINDRICAL
    dx = ymed(j)*cos(xmed(i))-Xplanet;
    dy = ymed(j)*sin(xmed(i))-Yplanet;
#ifdef Z
    dz = zmed(k)-Zplanet;
#endif
#endif
#ifdef SPHERICAL
    dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-Xplanet;
    dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-Yplanet;
#ifdef Z
    dz = ymed(j)*cos(zmed(k))-Zplanet;
#endif
#endif
  dist = sqrt(dx*dx + dy*dy + dz*dz);
  sink = exp(-pow(dist/rsink,4.));
  interm[ll] = dens[ll] * (gammasink*omega_bin*sink)*Vol(i,j,k);  // mdot[ll]; //*Vol(i,j,k); //dens[ll] * (gammasink*sink*sink)*Vol(i,j,k);
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}