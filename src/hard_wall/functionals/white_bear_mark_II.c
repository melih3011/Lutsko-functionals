#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

#define n3m 1.e-4

static double phi1(double n3) {

  if (n3<n3m) return 1.+n3*n3/9.;

  return 1.+(2.*n3-n3*n3+2.*(1.-n3)*log(1.-n3))/(3.*n3);

}

static double phi2(double n3) {
  
  if (n3<n3m) return 1.-4./9.*n3+n3*n3/18.;

  return 1.-(2.*n3-3.*n3*n3+2.*n3*n3*n3+2.*(1.-n3)*(1.-n3)*log(1.-n3))/(3.*n3*n3);

}

double pressure(double eta,double rho0,double a,double b)
{
  double hc;

  hc=1.0-eta;
  return rho0*(1+eta*(1+eta*(1-eta)))/(hc*hc*hc);
}

double sk0_bulk(double eta,double a,double b)
{
  double hc;

  hc=1.0-eta;
  return pow(hc,4)/(1.+4.*eta+4.*eta*eta-eta*eta*eta*(4.-eta));
}

double pressure_mixture(double n0,double n1,double n2,double n3)
{
  double hc;

  hc=(1.-n3)*(1.-n3)*(n3-1.);

  return (-(n2*n2*n2*(3.+(n3-2.)*n3))+12.*(n3-1.)*
	  (-3.*n0*(n3-1.)+n1*n2*(3 +n3*n3))*pi)/(36.*hc*pi);
}

double chemical_potential(double eta,double a,double b)
{
  double hc;

  hc=1.-eta;
  return eta*(8.+3.*(eta-3.)*eta)/(hc*hc*hc);
}

double detadmu(double eta, double a, double b) {

    return eta*sk0_bulk(eta, a, b);
}

double chemical_potential_mixture(double n0,double n1,double n2,double n3,
				  double hR)
{
  double hc,hc2,hc3;

  hc=1.-n3;
  hc2=hc*hc*(n3-1.);
  hc3=hR*n2-3*n3;

  return (hR*n3*(n2*(9.*(n3-5.)*hc*hc*n3*n3-9*hR*n2*(n3-1.)*n3*(1.+(n3-3.)*n3)+
		     hR*hR*n2*n2*(-2.+(n3-5.)*(n3-1.)*n3))- 
		 12.*hR*(n3-1.)*n3*(-3.*n1*(n3-5.)*(n3-1.)*n3+ 
				    hR*(n1*n2*(2.-6.*n3)+3.*n0*(n3-1.)*n3))*pi)
	  +(hR*n2-3.*n3)*hc2*(-3.*hR*n2*n3 + 9.*n3*n3 + 
			      2.*hR*hR*(n2*n2 - 12.*n1*n3*pi))*log(1.-n3))/
      (27.*hc2*n3*n3*n3);
}


double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2)
{
  return (-n0)*log(1-n3);
}

double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2)
{
    return phi1(n3)*(n1*n2-nv2*nv1)/(1.-n3);
}

double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b)
{
  double n13,hc;

  n13=1.-n3;
  hc=n13*n13;

  return phi2(n3)*n2*(n2*n2-3*nv2*nv2)/(24.*pi*n13*n13); 
}

double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR,double a,double b)
{
  double d21,d22,d23,ln3;

  ln3=1.-n3;
  
  d21=phi1(n3)*n1/ln3+phi2(n3)*(n2*n2-nv2*nv2)/(8.*pi*ln3*ln3);

  d22=phi1(n3)*n2/ln3;
      
  d23=(-log(ln3));

  return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);
}

double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR,double a,double b)
{
    double n32,n33,ln3,ln32,ln33,hr1,hr2;
    
    ln3=1.-n3;
    n32=n3*n3;
    n33 = pow(n3,3);
    ln32=ln3*ln3;
    ln33 = pow(1.-n3,3);
    hr1=n1*n2-nv1*nv2;
    hr2=n2*(n2*n2-3.*nv2*nv2);

    if (n3<n3m) return n0/ln3+hr1*(1.+20./9.*n3+7./2.*n32)+hr2/pi*(7./108.+13./72.*n3+7./20.*n32);

    return n0/ln3-hr1*2.*(n3-3.*n32+ln32*log(ln3))/(3.*ln32*n32)-
    hr2*(n3*(-2.+5.*n3-6.*n32+n33)-2.*ln33*log(ln3))/(36.*pi*ln33*n33);
}

double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	    double hR,double a,double b)
{
  double dv21,dv22,dv2,ln3;

  ln3=1.-n3;

  dv21=-phi1(n3)*nv1/ln3-phi2(n3)*n2*nv2/(4.*pi*ln3*ln3);

  dv22=-phi1(n3)*nv2/ln3;

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
  sprintf(fln,"%s","wb_II");
}
