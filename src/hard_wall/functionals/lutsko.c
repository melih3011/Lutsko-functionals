#include <math.h>
#include <stdio.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */
#define pi_by_2 pi/2

// Lutsko version 
// Two open parameters A and B must be included




// thermodynamic quantities
//=========================================
double z(double eta, double a, double b) {

	double z_py = (1.+eta+eta*eta)/pow(1.-eta, 3);

	return z_py + eta*eta/pow(1.-eta, 3)*((2./3.)*(4.*a+b)-3.);

}


double pressure(double eta, double rho0, double a, double b) {

	return z(eta, a, b)*rho0;
}


double chemical_potential(double eta, double a, double b) {

	return -log(1.-eta)+eta*(7.-11.*eta+4*eta*eta)/pow(1.-eta,3)+((4.*a+b)/3.)*(3.*eta*eta-eta*eta*eta)/pow(1.-eta,3);

}

double detadmu(double eta, double a, double b) {

	double num = eta*pow(1.-eta,4);
	double denom = 1.+4.*eta+(-5.+8.*a+2.*b)*eta*eta;

	return num/denom;
}


double surface_tension(double eta, double a, double b) {

	double hc = eta/(1.-eta);
	
	return (3./pi_4)*hc+3./(2.*pi_4)*pow(hc, 2)*(3.*a+2.*b);
}


double sk0_bulk(double eta, double a, double b) {

	double num = pow(1.-eta, 4);
	double denom = 1.+4.*eta+(2.*(4.*a+b)-5.)*eta*eta;

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
	double denom = 24.*pi*pow(1.-n3,2);

	return num/denom;
}
//=========================================





// myd-section
//=========================================

double myd2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {

	double d21, d22, d23, ln3=1.-n3;

	double trace1 = 2.*pow(nt1,2)+pow(nt2,2);
	double trace2 = 2.*pow(nt1,3)+pow(nt2,3);
	
	d21 = n1/ln3+((a+b)*pow(n2,2)-a*pow(n2v,2)-b*trace1)/(8.*pi*pow(ln3,2));
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
	d33 = ((a+b)*pow(n2,3)-3.*a*n2*pow(n2v,2)+3.*a*pow(n2v,2)*nt2-3.*b*n2*trace1+(2.*b-a)*trace2)/(12.*pi*pow(ln3,3));

	return d31+d32+d33;


}


double mydv2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double dv21, dv22, ln3=1.-n3;

	dv21 = -n1v/ln3+(a*n2v*(nt2-n2))/(pi_4*pow(ln3,2));
	dv22 = -n2v/ln3;

	return dv21+dv22/(pi_4*hR);

}

double mydt1(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double ln3 = 1.-n3;

	return nt1*((2.*b-a)*nt1-2.*b*n2)/(pi_4*pow(ln3,2));

}

double mydt2(double n0, double n1, double n2, double n3, double n1v, double n2v, double nt1, double nt2, double hR, double a, double b) {
	
	double ln3 = 1.-n3;

	return (a*pow(n2v,2)-2.*b*n2*nt2+(2.*b-a)*pow(nt2,2))/(8.*pi*pow(ln3,2));

}

//=========================================




void filename(char* fln) {
	sprintf(fln, "%s", "lk");
}
