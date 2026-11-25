#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

double pressure(double eta,double rho0,double a,double b)
{
  double hc;

  hc=1.0-eta;
  return rho0*(1+eta+eta*eta)/(hc*hc*hc);
}

double chemical_potential(double eta, double a,double b)
{
  double hc;

  hc=eta/(1-eta);
  return hc*(7.0+hc*(7.5+3*hc))-log(1-eta);
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
  double lz,hc,ln3;

  lz= (n2==0) ? 0 : nv2/n2;
  hc=n2*(1-lz*lz);
  ln3=1.-n3;
  return hc*hc*hc/(24.0*pi*ln3*ln3);
}

double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,
	    double hR)
{
  double d21,d22,d23,ln3,lz,hc;

  lz= (n2==0) ? 0 : nv2/n2;
  hc=n2*(1-lz*lz);
  ln3=1.-n3;

  d21=n1/ln3+hc*hc*(1+lz*lz)/(8*pi*ln3*ln3);

  d22=n2/ln3;

  d23=(-log(ln3));

  return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);
}

double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,
	    double hR)
{
  double ln3,lz,hc;

  lz= (n2==0) ? 0 : nv2/n2;
  ln3=1.-n3;
  hc=n2*(1-lz*lz)/ln3;

  return n0/ln3+(n1*n2-nv1*nv2)/(ln3*ln3)+hc*hc*hc/(12*pi);
}

double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,
	    double hR)
{
  double dv21,dv22,dv2,ln3,lz,hc;

  lz= (n2==0) ? 0 : nv2/n2;
  ln3=1.-n3;
  hc=(1-lz*lz);

  dv21=nv1/ln3+n2*nv2*hc*hc/(pi_4*ln3*ln3);

  dv22=nv2/ln3;

  return (dv21+dv22/(pi_4*hR));
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
  sprintf(fln,"%s","rfq3");
}
