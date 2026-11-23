#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
/* Test-only analytic reference and Fourier sums; no Python numerical kernel. */
int main(int argc,char **argv) {
 if(argc!=3 && argc!=4 && argc!=5)return 2;
 double background=argc==5?strtod(argv[4],NULL):1.;
 if(!isfinite(background)||background<=0)return 2;
 FILE *fp=fopen(argv[1],"r");if(!fp)return 3;
 double t,reference,sample,d=20,c=137.035999084;
 const double frequencies[3]={6,10,14};
 double complex ref[3]={0},mat[3]={0};
 int rows=0;
 while(fscanf(fp,"%lf %lf %lf",&t,&reference,&sample)==3) {
  for(int k=0;k<3;k++) {double complex phase=cexp(I*frequencies[k]*t);ref[k]+=reference*phase;mat[k]+=sample*phase;}
  rows++;
 }
 fclose(fp);if(rows<100)return 4;
 fp=fopen(argv[2],"w");if(!fp)return 5;
 double worst=0;
 for(int k=0;k<3;k++) {
  double w=frequencies[k];
  double complex epsilon=argc>=4 ? background+400./(625.-w*w-2.*I*w) : 2.25;
  double complex n=csqrt(epsilon),r=(n-1)/(n+1);
  double complex exact=4*n/((1+n)*(1+n))*cexp(I*(n-1)*w*d/c)/(1-r*r*cexp(2*I*n*w*d/c));
  double complex measured=mat[k]/ref[k];double error=cabs(measured-exact)/cabs(exact);
  if(!isfinite(error)||cabs(ref[k])<1e-12)return 6;
  if(error>worst)worst=error;
  fprintf(fp,"%.17g %.17g %.17g %.17g %.17g %.17g\n",w,creal(measured),cimag(measured),creal(exact),cimag(exact),error);
 }
 fclose(fp);printf("maximum complex transmission relative error %.17g\n",worst);
 return worst>.02?7:0;
}
