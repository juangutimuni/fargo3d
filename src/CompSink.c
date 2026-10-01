//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void ComputeSink (real dt){

int nb, i, j;
  nb = Sys->nb;
  for (i = 0; i < nb; i++) {
    if (ThereArePlanets) {
      Xplanet = Sys->x[i];
      Yplanet  = Sys->y[i];
      Zplanet  = Sys->z[i];
      VXplanet = Sys->vx[i];
      VYplanet = Sys->vy[i];
      VZplanet = Sys->vz[i];
      MplanetVirtual = Sys->mass[i];
    } else {
      Xplanet = Yplanet = Zplanet = VXplanet = VYplanet = VZplanet = MplanetVirtual = 0.0;
    } 
    FARGO_SAFE(_ComputeSink(dt));
    }

}
void _ComputeSink_cpu (real dt) {

//<USER_DEFINED>
  INPUT(Density);
  OUTPUT(Density);
  INPUT(Vx_temp);
  INPUT(Vy_temp);
  OUTPUT(Vx_temp);
  OUTPUT(Vy_temp);
  #ifdef Z
  INPUT(Vz_temp);
  OUTPUT(Vz_temp);
  #endif
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* vx = Vx_temp->field_cpu;
  real* vy = Vy_temp->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ;
  real gammasink = GAMMASINK;
  real rsink = RSINK;
  real dsink = DELTASINK;
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
  real phip;
  real alpha;
  real vpp;
  real vpr;
  real vrstar;
  real vpstar;
  real vrcstar;
  real vpcstar;
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
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
  #ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-Xplanet;
	dy = ymed(j)*sin(xmed(i))-Yplanet;
#ifdef Z
	dz = zmed(k)-Zplanet;
#endif
  phip = atan2(Yplanet, Xplanet);
  alpha = phip - xmed(i);
  vpr = VXplanet * cos(phip) + VYplanet * sin(phip); 
  vpp = - VXplanet * sin(phip) + VYplanet * cos(phip);    
  vpcstar = vx[ll]; //0.25*(vx[ll]+vx[lym]+vx[lxp]+vx[lxp-pitch]); 
  vrcstar = vy[ll]; //0.25*(vy[ll]+vy[lyp]+vy[lxm]+vy[lxm+pitch]);
  vrstar = vpr + (vy[ll] - vpr)*(cos(alpha)*cos(alpha) + dsink *sin(alpha)*sin(alpha)+(vpcstar-vpp)*(1-dsink)*sin(alpha)*cos(alpha)); 
  vpstar = vpp + (vrcstar - vpr)*(1-dsink)*sin(alpha)*cos(alpha) + (vx[ll] - vpp)*(dsink*cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha));
  #endif
#ifdef SPHERICAL
	dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-Xplanet;
	dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-Yplanet;
#ifdef Z
	dz = ymed(j)*cos(zmed(k))-Zplanet;
#endif
#endif
  dist = sqrt(dx*dx + dy*dy + dz*dz);
  sink = 1 - dist*dist /(rsink*rsink);
	if(dist/rsink < 1.){
    vx[ll] = vx[ll] -gammasink*0.5*(dens[ll]+dens[lxm])*dt*sink*sink*vpstar;
    vy[ll] = vy[ll] -gammasink*0.5*(dens[ll]+dens[lym])*dt*sink*sink*vrstar;
#ifdef Z
    vz[ll] = vz[ll] * (1.-gammasink*0.5*(dens[ll]+dens[lzm])*dt*sink*sink);
#endif
    dens[ll] = dens[ll] * (1.-gammasink*dt*sink*sink);
  } 
    
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
