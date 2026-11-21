#include <math.h>
#include <stddef.h>
/* Validation-only reductions/collocation; production solver is unchanged. */
int energy_center_h(int nt,int nx,const double *raw,double *out){
 if(nt<1||nx<2)return 1;
 for(int t=0;t<nt;t++)for(int x=1;x<nx;x++)for(int k=0;k<6;k++){
  size_t p=((size_t)t*nx+x)*6+k,q=((size_t)t*(nx-1)+x-1)*6+k;
  if(!isfinite(raw[p])||!isfinite(raw[p-6]))return 2;
  out[q]=(k==4||k==5)?0.5*(raw[p]+raw[p-6]):raw[p];
 }return 0;
}
int energy_split_x(int n,const double *fields,double *ep,double *hp,double *em,double *hm){
 if(n<1)return 1;
 for(int t=0;t<n;t++){
  const double *f=fields+6*t;double *a=ep+3*t,*b=em+3*t,*h=hp+3*t,*g=hm+3*t;
  for(int k=0;k<6;k++)if(!isfinite(f[k]))return 2;
  a[0]=b[0]=h[0]=g[0]=0;
  a[1]=.5*(f[1]+f[5]);a[2]=.5*(f[2]-f[4]);
  b[1]=.5*(f[1]-f[5]);b[2]=.5*(f[2]+f[4]);
  h[1]=-a[2];h[2]=a[1];g[1]=b[2];g[2]=-b[1];
 }return 0;
}
int energy_field_sum(int n,const double *fields,double volume,double *out){
 if(n<1||!isfinite(volume)||volume<=0)return 1;
 double v=0;
 for(int i=0;i<n;i++)for(int k=0;k<6;k++){
  double a=fields[6*i+k];if(!isfinite(a))return 2;v+=a*a;
 }
 *out=v*volume/(8*acos(-1.));return isfinite(*out)?0:2;
}
int energy_weighted_sum(int n,const double *values,const double *weights,double *out){
 if(n<1)return 1;double s=0;
 for(int i=0;i<n;i++){
  if(!isfinite(values[i])||!isfinite(weights[i]))return 2;s+=values[i]*weights[i];
 }*out=s;return isfinite(s)?0:2;
}

/* Four x planes: left-neighbor,left,right-neighbor,right. z invariant. */
int energy_periodic_planes(int nt,int ny,const double *raw,double *out){
 if(nt<1||ny<2)return 1;
 for(int t=0;t<nt;t++)for(int side=0;side<2;side++)for(int y=0;y<ny;y++){
  int plane=2*side+1,ym=(y+ny-1)%ny;
  const double *a=raw+(((size_t)t*4+plane)*ny+y)*6;
  const double *b=raw+(((size_t)t*4+plane-1)*ny+y)*6;
  const double *ay=raw+(((size_t)t*4+plane)*ny+ym)*6;
  const double *by=raw+(((size_t)t*4+plane-1)*ny+ym)*6;
  double *v=out+(((size_t)t*2+side)*ny+y)*6;
  for(int k=0;k<6;k++)if(!isfinite(a[k])||!isfinite(b[k])||!isfinite(ay[k])||!isfinite(by[k]))return 2;
  v[0]=.5*(a[0]+b[0]);v[1]=.5*(a[1]+ay[1]);v[2]=a[2];
  v[3]=.5*(a[3]+ay[3]);v[4]=.5*(a[4]+b[4]);v[5]=.25*(a[5]+b[5]+ay[5]+by[5]);
 }return 0;
}
/* Fixed test profile: nonmagnetic, isotropic dielectric, TE, periodic y/z.
   Recover H at the E time from the exact Yee H update (region is outside PML).
   raw arrays use C xyz order; epsilon is the z-edge dielectric background. */
int energy_yee_te_end(int nx,int ny,int nz,int first,int last,double dx,double dy,double dt,
                     double volume,const double *ez,const double *hx,const double *hy,
                     const double *epsilon,double *stats){
 if(nx<3||ny<2||nz<1||first<1||last>=nx-1||last<first||dx<=0||dy<=0||dt<=0||volume<=0)return 1;
 double energy=0,raw_energy=0;
 for(int x=first;x<=last;x++)for(int y=0;y<ny;y++)for(int z=0;z<nz;z++){
  size_t i=((size_t)x*ny+y)*nz+z,ix=i+(size_t)ny*nz,iy=((size_t)x*ny+(y+1)%ny)*nz+z;
  if(!isfinite(ez[i])||!isfinite(ez[ix])||!isfinite(ez[iy])||!isfinite(hx[i])||!isfinite(hy[i])||
     !isfinite(epsilon[i])||epsilon[i]<=0)return 2;
  double hxc=hx[i]+.5*137.035999084*dt*(ez[iy]-ez[i])/dy;
  double hyc=hy[i]-.5*137.035999084*dt*(ez[ix]-ez[i])/dx;
  double w=(x==first||x==last)?.5:1.;
  energy+=w*(epsilon[i]*ez[i]*ez[i]+hxc*hxc+hyc*hyc);
  raw_energy+=w*(epsilon[i]*ez[i]*ez[i]+hx[i]*hx[i]+hy[i]*hy[i]);
 }
 stats[0]=energy*volume/(8*acos(-1.));stats[1]=raw_energy*volume/(8*acos(-1.));
 return isfinite(stats[0])&&isfinite(stats[1])?0:2;
}
