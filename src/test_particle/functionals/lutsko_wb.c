#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

#define min_eta 0.01

// Lutsko version 
// Two open parameters A and B must be included


//f_wb defining the prefactor of phi3 in white-bear version
//=========================================
double f_wb(double eta) {

	double num = eta + pow(1.-eta,2)*log(1.-eta);
	double denom = 36.*pi*pow(eta*(1.-eta),2);

	if (eta <= min_eta) return 1./(24.*pi)+2.*eta/(27.*pi)+5.*eta*eta/(48.*pi);

	return num/denom;

}

double df_wb(double eta) {

	double num = eta*(2.-5*eta+eta*eta) + 2.*pow(1.-eta,3)*log(1.-eta);
	double denom = 36.*pi*pow(eta*(1.-eta),3);

	if (eta <= min_eta) return 2./(27.*pi)+5.*eta/(24.*pi)+2.*eta*eta/(5.*pi);

	return -num/denom;

}
//=========================================


// thermodynamic quantities
//=========================================
//compressibility is derived from dphi/dn3|b
double z(double eta, double a, double b) {

	return 1./(1.-eta) + 3.*eta/pow(1.-eta,2) + df_wb(eta)*8.*pi*eta*eta*(4.*a+b);
}

//pressure derived from thermodynamic relation bp = -bf + b eta df/deta
double pressure(double eta, double rho0, double a, double b) {

	return ((1.+eta-2.*eta*eta)/pow(1.-eta,3)+(2./3.)*(4.*a+b)*eta*eta/pow(1.-eta,3)*(1.-eta/3.))*rho0;

}


/*
double surface_tension(double eta, double a, double b) {

	double hc = eta/(1.-eta);
	
	return (3./pi_4)*hc + f_wb(eta)*9.*eta*eta*(3.*a+2.*b);
}
*/

//chemical potential derived from thermodynamic relation bmu = df/drho
double chemical_potential(double eta, double a, double b) {

	double mu = -eta/(9.*pow(1.-eta,3))*(2.*(4.*a+b)*(-1.-2.*eta+eta*eta)-9.*(7.-11.*eta+4.*eta*eta))+(-9.+8.*a+2.*b)*log(1.-eta)/9.;

	return mu;
}

double detadmu(double eta, double a, double b) {

	double num = 9.*eta*pow(1.-eta,4);
	double denom = (9.+36.*eta+9.*(-5.+8.*a+2.*b)*eta*eta+2.*(4.*a+b)*pow(eta,3)*(eta-4.));

	return num/denom;

}

double sk0_bulk(double eta, double a, double b) {

	double num = 9.*pow(1.-eta, 4);
	double denom = 9.+36.*eta+9.*(-5.+8.*a+2.*b)*eta*eta+2.*(4.*a+b)*pow(eta,3)*(eta-4.);

	return num/denom;

}
//=========================================



// myphi-section
//=========================================

double myphi1(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2) {

	return -n0*log(1.-n3);
}

double myphi2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2) {

	return (n1*n2-n2v*n1v)/(1.-n3);
}

double myphi3(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double a, double b) {

	double trace1 = 2.*pow(nt1,2)+pow(nt2,2);
	double trace2 = 2.*pow(nt1,3)+pow(nt2,3);
	
	double num = (a+b)*pow(n2,3)-3.*a*n2*pow(n2v,2)+3.*a*pow(n2v,2)*nt2-3.*b*n2*trace1+(2.*b-a)*trace2;

	return f_wb(n3)*num;
}
//=========================================





// myd-section
//=========================================

double myd2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {

	double d21, d22, d23, ln3=1.-n3;

	double trace1 = 2.*pow(nt1,2)+pow(nt2,2);
	double trace2 = 2.*pow(nt1,3)+pow(nt2,3);
	
	d21 = n1/ln3+f_wb(n3)*3.*((a+b)*pow(n2,2)-a*pow(n2v,2)-b*trace1);
	d22 = n2/ln3;
	d23 = -log(ln3);
	
	return d21+d22/(pi_4*hR)+d23/(pi_4*hR*hR);

}

double myd3(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double d31,d32,d33,ln3=1.-n3;
	
	double trace1 = 2.*pow(nt1,2)+pow(nt2,2);
	double trace2 = 2.*pow(nt1,3)+pow(nt2,3);

	d31 = n0/ln3;
	d32 = (n1*n2-n1v*n2v)/pow(ln3,2);
	d33 = ((a+b)*pow(n2,3)-3.*a*n2*pow(n2v,2)+3.*a*pow(n2v,2)*nt2-3.*b*n2*trace1+(2.*b-a)*trace2)*df_wb(n3);

	return d31+d32+d33;


}


double mydv2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double dv21, dv22, ln3=1.-n3;

	dv21 = -n1v/ln3+6.*a*n2v*(nt2-n2)*f_wb(n3);
	dv22 = -n2v/ln3;

	return dv21+dv22/(pi_4*hR);

}

double mydt1(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double ln3 = 1.-n3;

	//return 3.*nt1*((2.*b-a)*nt1-4.*b*n2)*f_wb(n3);
	return 6.*f_wb(n3)*nt1*((2.*b-a)*nt1-2.*b*n2);

}

double mydt2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double ln3 = 1.-n3;

	return 3.*(a*pow(n2v,2)-2.*b*n2*nt2+(2.*b-a)*pow(nt2,2))*f_wb(n3);
}

//=========================================




void filename(char* fln) {
	sprintf(fln, "%s", "lkwb");
}
