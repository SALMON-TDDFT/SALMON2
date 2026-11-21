#define _XOPEN_SOURCE 700
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Independent test-only outgoing-wave cylinder reference; no production solver. */
static double complex cj0(double complex z) {
 double complex term=1.,sum=1.;
 for(int k=1;k<100;k++){term*=-(z*z/4.)/(k*(double)k);sum+=term;if(cabs(term)<1e-17*cabs(sum))break;}
 return sum;
}
static double complex cj1(double complex z) {
 double complex term=z/2.,sum=term;
 for(int k=1;k<100;k++){term*=-(z*z/4.)/(k*(double)(k+1));sum+=term;if(cabs(term)<1e-17*cabs(sum))break;}
 return sum;
}
int main(int argc,char **argv) {
 if(argc!=8)return 2;
 double w=strtod(argv[3],NULL),radius=strtod(argv[4],NULL),wp=strtod(argv[5],NULL);
 double gamma=strtod(argv[6],NULL),amplitude=strtod(argv[7],NULL);
 if(!isfinite(w)||!isfinite(radius)||!isfinite(wp)||!isfinite(gamma)||!isfinite(amplitude))return 2;
 if(w<=0||radius<=0||wp<=0||gamma<=0)return 2;
 FILE *fp=fopen(argv[1],"r");if(!fp)return 3;
 double t,e,old_t=0,old_e=0;double complex integral=0;int rows=0;
 while(fscanf(fp,"%lf %lf",&t,&e)==2){
  if(!isfinite(t)||!isfinite(e)||t<old_t){fclose(fp);return 4;}
  integral+=(t-old_t)*.5*(old_e*cexp(I*w*old_t)+e*cexp(I*w*t));old_t=t;old_e=e;rows++;
 }
 fclose(fp);if(rows<100)return 4;
 double complex denominator=-w*w-I*gamma*w,eps=1.+wp*wp/denominator,n=csqrt(eps);
 double x=w*radius/137.035999084;
 double complex h0=j0(x)+I*y0(x),hprime=-j1(x)-I*y1(x);
 double complex analytic=wp*wp*amplitude/(denominator+wp*wp) *
   (1.+hprime/(-n*cj1(n*x)*h0-cj0(n*x)*hprime));
 double error=cabs(integral-analytic)/cabs(analytic);
 if(!isfinite(error))return 5;
 fp=fopen(argv[2],"w");if(!fp)return 3;
 fprintf(fp,"{\"omega\":%.17g,\"measured\":[%.17g,%.17g],\"analytic\":[%.17g,%.17g],\"relative_error\":%.17g,\"duration\":%.17g}\n",
   w,creal(integral),cimag(integral),creal(analytic),cimag(analytic),error,old_t);
 fclose(fp);printf("cylinder-center E_z complex relative error %.17g\n",error);return 0;
}
