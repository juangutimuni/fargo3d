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
  real* vz = Vz_temp->field_cpu;
  real* mdot = Mdot->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY-1;
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
  real phip;
  real alpha;
  real vpp, vpt, vpr;
  real vrstar, vpstar, vtstar;
  real vrcstar, vpcstar, vtcstar;
  real vrc, vpc, vtc;
  real thetap;
  real ar, ap, at;
  real wr, wp, wt;
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
    for (j=1; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;
  mdot[ll] = 0;
#ifdef CYLINDRICAL
	dx = ymed(j)*cos(xmed(i))-Xplanet;
	dy = ymed(j)*sin(xmed(i))-Yplanet;
#ifdef Z
	dz = zmed(k)-Zplanet;
#endif
  phip = atan2(ymed(j)*sin(xmed(i))-Yplanet, ymed(j)*cos(xmed(i))-Xplanet);   
  alpha = phip - xmed(i);
  vpr =  VXplanet * cos(xmed(i)) + VYplanet * sin(xmed(i)); 
  vpp = -VXplanet * sin(xmed(i)) + VYplanet * cos(xmed(i));
  vrc = vy[ll];
  vpc = vx[ll];  
  vpcstar = 0.25*(vx[ll]+vx[lym]+vx[lxp]+vx[lxp-pitch]); 
  vrcstar = 0.25*(vy[ll]+vy[lyp]+vy[lxm]+vy[lxm+pitch]);
  vrstar = vpr + (vy[ll] - vpr)*(cos(alpha)*cos(alpha) + dsink *sin(alpha)*sin(alpha))+(vpcstar-vpp)*(1.-dsink)*sin(alpha)*cos(alpha); 
  vpstar = vpp + (vrcstar - vpr)*(1.-dsink)*sin(alpha)*cos(alpha) + (vx[ll] - vpp)*(dsink*cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha));
  #endif
#ifdef SPHERICAL
	dx = ymed(j)*cos(xmed(i))*sin(zmed(k))-Xplanet;
	dy = ymed(j)*sin(xmed(i))*sin(zmed(k))-Yplanet;
#ifdef Z
	dz = ymed(j)*cos(zmed(k))-Zplanet;
  phip = atan2(dy, dx);   
  thetap = acos(fmax(-1., fmin(1., dz/(sqrt(dx*dx + dy*dy + dz*dz) + 1e-30))));
  alpha = phip - xmed(i);
  vpr =  VXplanet * cos(xmed(i)) * sin(zmed(k)) + VYplanet * sin(xmed(i)) * sin(zmed(k)) + VZplanet * cos(zmed(k));
  vpp = -VXplanet*sin(xmed(i)) + VYplanet*cos(xmed(i));
  vpt =  VXplanet * cos(xmed(i)) * cos(zmed(k)) + VYplanet * sin(xmed(i)) * cos(zmed(k)) - VZplanet * sin(zmed(k));
  ar = sin(zmed(k))*sin(thetap)*cos(alpha) + cos(zmed(k))*cos(thetap);
  ap = sin(thetap)*sin(alpha);
  at = cos(zmed(k))*sin(thetap)*cos(alpha) - sin(zmed(k))*cos(thetap);
  vrc = 0.5*(vy[ll] + vy[lyp]);   
  vpc = 0.5*(vx[ll] + vx[lxp]);   
  vtc = 0.5*(vz[ll] + vz[lzp]);   
  wr = vrc - vpr;
  wp = vpc - vpp;
  wt = vtc - vpt;
  vrstar = vpr + (dsink + (1-dsink)*ar*ar)*wr+ (1-dsink)*ar*ap*wp + (1-dsink)*ar*at*wt;
  vpstar = vpp + (1-dsink)*ap*ar*wr + (dsink + (1-dsink)*ap*ap)*wp + (1-dsink)*ap*at*wt; 
  vtstar = vpt + (1-dsink)*at*ar*wr + (1-dsink)*at*ap*wp + (dsink + (1-dsink)*at*at)*wt;
#endif
#endif
  dist = sqrt(dx*dx + dy*dy + dz*dz);
  sink = exp(-pow(dist/rsink,4.));
  // the staggering needs to be checked
  vx[ll] = vx[ll] - gammasink*omega_bin*dt*sink*(vpstar - vpc); // vx[ll] -gammasink*0.5*(dens[ll]+dens[lxm])*dt*omega_bin*sink*vpstar;
  vy[ll] = vy[ll] - gammasink*omega_bin*dt*sink*(vrstar - vrc); // vy[ll] -gammasink*0.5*(dens[ll]+dens[lym])*dt*omega_bin*sink*vrstar;
#ifdef Z
  vz[ll] = vz[ll] - gammasink*omega_bin*dt*sink*(vtstar - vtc); // vz[ll] * (1.-gammasink*0.5*(dens[ll]+dens[lzm])*dt*omega_bin*sink);
#endif
  dens[ll] = dens[ll] * (1.-gammasink*dt*omega_bin*sink); // mdot[ll] = dens[ll] * (gammasink*omega_bin*sink)*Vol(i,j,k);
    
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
