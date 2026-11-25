#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

double pressure(double eta,double rho0, double a, double b)
{
  double hc;

  hc=1.0-eta;
  return rho0*(1+eta*(1+eta*(1-eta)))/(hc*hc*hc);
}

double pressure_mixture(double n0,double n1,double n2,double n3)
{
  double hc;

  hc=1/(1.-n3);

  return hc*(n0+hc*(n1*n2+hc*(n2*n2*n2*(1.-n3/3.))/(12.*pi)));
}

double chemical_potential(double eta, double a, double b)
{
  double hc;

  hc=1.-eta;
  return eta*(8.+3.*(eta-3.)*eta)/(hc*hc*hc);
}

double chemical_potential_mixture(double n0,double n1,double n2,double n3,
				  double hR)
{
  double hc;

  hc=1.-n3;

  return -log(hc)+hR*n2/hc+pi_4*hR*hR*(n1/hc+n2*n2*(n3+hc*hc*log(hc))/
				       (12.*pi*hc*hc*n3*n3))+
    pi_4*hR*hR*hR/3.0*(n0/hc+n1*n2/(hc*hc)-
		       n2*n2*n2*(n3*(2.-5.*n3+n3*n3)+2*hc*hc*hc*log(hc))
		       /(36.*pi*hc*hc*hc*n3*n3*n3));
}

double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2)
{
  return (-n0)*log(1-n3);
}

double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2)
{
    return (n1*n2-nv2*nv1)/(1-n3);
}

double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b)
{
  double n13,hc,hc2,lz;

  if (n3==0) return 0;
  lz= (n2==0) ? 0 : nv2/n2;

  n13=1.-n3;
  hc=n13*n13;
  hc2=n2*(1-lz*lz);

  return hc2*hc2*hc2*(n3+hc*log(n13))/(36.0*pi*n3*n3*hc); 
}

double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR, double a,double b)
{
  double d21,d22,d23,hc,hc2,ln3,lz;

  if (n3==0.) return 0;
  lz= (n2==0) ? 0 : nv2/n2;

  ln3=1.-n3;

  hc=ln3*n3;
  hc2=n2*(1.-lz*lz);

  d21=n1/ln3+hc2*hc2*(1+lz*lz)*(n3+ln3*ln3*log(ln3))/(12.*pi*hc*hc);

  d22=n2/ln3;

  d23=(-log(ln3));

  return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);
}

double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR, double a,double b)
{
  double o1,o2,o3,o4,o5,o6,o7,lz;

  if (n3==0.) return 0;

  lz= (n2==0) ? 0 : nv2/n2;

  o1=-1.+n3;
  o2=o1*o1*o1;
  o3=n3*n3*n3;
  o4=n2*(1.-lz*lz);
  o7=o4*o4*o4;
  
  return (-36*pi*o1*(-(n1*n2)+nv1*nv2+n0*o1)*o3+
	  n3*(2.+(-5.+n3)*n3)*o7-
	  2.*o2*o7*log(1.-n3))/(36*pi*o2*o3);
}

double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR, double a,double b)
{
  double dv21,dv22,dv2,ln3,hc,hc2,lz;

  if (n3==0.) return 0;

  lz= (n2==0) ? 0 : nv2/n2;
  ln3=1.-n3;
  
  hc=ln3*n3;
  hc2=(1.-lz*lz);

  dv21=nv1/ln3+n2*nv2*hc2*hc2*(n3+ln3*ln3*log(ln3))/(6.*pi*hc*hc);

  dv22=nv2/ln3;

  return -(dv21+dv22/(pi_4*hR));
}
double mydt1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
      double hR, double a, double b) {

        return 0.;
}


double mydt2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
      double hR, double a, double b) {

        return 0.;
}
void filename(char *fln)
{
  sprintf(fln,"%s","wbq3");
}
