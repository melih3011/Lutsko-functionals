#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769
#define pi_4 12.5663706143591729539

#define N 4096
#define R 100

#define ITS 40000


double rho[N], n3[N], n2[N], n2v[N], nt1[N], nt2[N];
double w3[2*R+1][2*R+1], w2[2*R+1][2*R+1], w2v[2*R+1][2*R+1];
double phi[N],phiid[N],dphi[N],dep[N],ads[N],Vext[N];

double eta,rho0,mu,xt,dz=1.0/(1.0*R),a1,p_bulk,dedmu,dgdeta,Gamma;
double A,B;

FILE *fp;

double get_magnitude(const double* arr) {

  int i;
  double max=fabs(arr[0]);

  for (i=1; i<N-2*R; ++i) {

    if (fabs(arr[i]) > max) max = fabs(arr[i]);
  }

  return max;
  
}

double integrate_array_section(const double *f, double dx, int i0, int i1)
{
  double sum = 0.0;
  int  i;
    
  if (i0 == i1) return 0.0;
  else
    {
      // Eq.(4.1.14) of Numerical Recipes
      // assume n>6
      
      sum=3./8.*(f[i0]+f[i1])+7./6.*(f[i0+1]+f[i1-1])+23./24*(f[i0+2]+f[i1-2]);

      for (i=i0+3;i<=i1-3;i++) sum+=f[i];

      return sum*dx;
    }
}  


void param_init() {

  double pressure(double eta,double rho0,double a,double b);
  double chemical_potential(double eta,double a,double b);
  double detadmu(double eta, double a, double b);
  double sk0_bulk(double eta,double a,double b);
  
  rho0=3.0*eta/(pi_4);
  mu=chemical_potential(eta,A,B);
  xt=sk0_bulk(eta,A,B);
  p_bulk=pressure(eta,rho0,A,B);
  dedmu=detadmu(eta,A,B);

}

void init() {
	double lw3[2*R+1], lw2[2*R+1], lw2v[2*R+1];
	double z,lR=R*dz;

	param_init();

  for (int i=0; i<=2*R; i++) {
    z = (i-R)*dz;
    lw3[i] = pi*dz*(lR*lR-z*z);
    lw2[i] = pi_2*dz*lR;
    lw2v[i] = pi_2*dz*z;
  }
	
		
	// initialize weights in set arrays
	// these weights are obtained from textbook on numerical integration
	double p1[2], p2[3], p3[4], p4[5];
	p1[0] = 0.5; p1[1] = 0.5;
	p2[0] = 1./3.; p2[1] = 4./3.; p2[2] = 1./3.;
	p3[0] = 3./8.; p3[1] = 9./8; p3[2] = 9./8.; p3[3] = 3./8.;
	p4[0] = 14./45.; p4[1] = 64./45.; p4[2] = 24./45.; p4[3] = 64./45.; p4[4] = 14./45.;


	for (int i=0; i<5; i++) {
		for (int j=0; j<2*R; j++) {
			w3[i][j] = 0.0; w2[i][j] = 0.0; w2v[i][j] = 0.0;
		}
	}
	
  for (int j=0; j<2; j++) {
    w3[1][j] = lw3[j]*p1[j]; w2[1][j] = lw2[j]*p1[j]; w2v[1][j] = lw2v[j]*p1[j];
	}

	for (int j=0; j<3; j++) {
		w3[2][j] = lw3[j]*p2[j]; w2[2][j] = lw2[j]*p2[j]; w2v[2][j] = lw2v[j]*p2[j];
	}

	for (int j=0; j<4; j++) {
		w3[3][j] = lw3[j]*p3[j]; w2[3][j] = lw2[j]*p3[j]; w2v[3][j] = lw2v[j]*p3[j];
	}

	for (int j=0; j<5; j++) {
		w3[4][j] = lw3[j]*p4[j]; w2[4][j] = lw2[j]*p4[j]; w2v[4][j] = lw2v[j]*p4[j];
	}
	
	for (int i=5; i<=2*R; i++) {
		for (int j=0; j<2*R-i; j++) {
			w3[i][2*R-j] = 0.0; w2[i][2*R-j] = 0.0; w2v[i][2*R-j] = 0.0;
		}
	
	w3[i][i] = 3./8.*lw3[i]; w3[i][i-1] = 7./6.*lw3[i-1]; w3[i][i-2] = 23./24.*lw3[i-2]; 
	w2[i][i] = 3./8.*lw2[i]; w2[i][i-1] = 7./6.*lw2[i-1]; w2[i][i-2] = 23./24.*lw2[i-2]; 
	w2v[i][i] = 3./8.*lw2v[i]; w2v[i][i-1] = 7./6.*lw2v[i-1]; w2v[i][i-2] = 23./24.*lw2v[i-2]; 

	for (int j=3; j<=i-3; j++) {
		w3[i][j] = lw3[j]; w2[i][j] = lw2[j]; w2v[i][j] = lw2v[j];
	} 

	w3[i][2] = 23./24.*lw3[2]; w3[i][1] = 7./6.*lw3[1]; w3[i][0] = 3./8.*lw3[0];
	w2[i][2] = 23./24.*lw2[2]; w2[i][1] = 7./6.*lw2[1]; w2[i][0] = 3./8.*lw2[0];
	w2v[i][2] = 23./24.*lw2v[2]; w2v[i][1] = 7./6.*lw2v[1]; w2v[i][0] = 3./8.*lw2v[0];
	}

}


void weighted_densities() {

	int i, j, jp, jm;
	double n3i, n2i, n2vi;

	for (i=0; i<N; i++) {
		n3[i] = 0.0;
		n2[i] = 0.0;
		n2v[i] = 0.0;
	}

	for (i=0; i<N-R; ++i) {

		n3i = 0.0;
		n2i = 0.0;
		n2vi = 0.0;

		jp = i+R;
		jm = i-R;

		if (jm<R) jm = R;

		for (j=jm; j<=jp; ++j) {

			n3i += rho[j]*w3[jp-jm][R+i-j];
			n2i += rho[j]*w2[jp-jm][R+i-j];
			n2vi += rho[j]*w2v[jp-jm][R+i-j];
		}

		n3[i] = n3i;
		n2[i] = n2i;
		n2v[i] = n2vi;
		nt1[i] = n3i;
		nt2[i] = n2i-2.*n3i;
	}

	n3[0]=0.0; n2[0]=0.0; n2v[0]=0.0;

}


void derivative() {
  double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	      double hR,double a,double b);
  double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	      double hR,double a,double b);
  double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	       double hR,double a,double b);
  double mydt1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	       double hR,double a,double b);
  double mydt2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,
	       double hR,double a,double b);  

  int i,j,jp,jm;
  double d2[N],d3[N],dv2[N],dt1[N],dt2[N];
  double lR=R*dz;
  double c3i,c2i,c2vi,ct1,ct2;
  double ln3,ln2,ln1,ln0,ln2v,ln1v,lnt1,lnt2;

  for (i=0;i<N-R;i++) {
    ln3=n3[i];
    ln2=n2[i];      
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    ln2v=n2v[i];
    ln1v=ln2v/(pi_4*lR);
    lnt1=nt1[i];
    lnt2=nt2[i];  

    d2[i]=myd2(ln0,ln1,ln2,ln3,ln1v,ln2v,lnt1,lnt2,lR,A,B);
    d3[i]=myd3(ln0,ln1,ln2,ln3,ln1v,ln2v,lnt1,lnt2,lR,A,B);
    dv2[i]=mydv2(ln0,ln1,ln2,ln3,ln1v,ln2v,lnt1,lnt2,lR,A,B);
    dt1[i]=mydt1(ln0,ln1,ln2,ln3,ln1v,ln2v,lnt1,lnt2,lR,A,B);
    dt2[i]=mydt2(ln0,ln1,ln2,ln3,ln1v,ln2v,lnt1,lnt2,lR,A,B);

  }

  for (i=N-R;i<N;i++) {
    d3[i]=0.0;
    d2[i]=0.0;
    dv2[i]=0.0;
    dt1[i]=0.0;
    dt2[i]=0.0;
  }

	for (i=0; i<N-R; ++i) {

		c3i = 0.0;
		c2i = 0.0;
		c2vi = 0.0;
		ct1 = 0.0;
		ct2 = 0.0;

		jp = i+R;
		jm = i-R;

		if (jm<0) jm=0;

		for (j=jm; j<=jp; ++j) {

			c3i += d3[j]*w3[jp-jm][R+i-j];
			c2i += d2[j]*w2[jp-jm][R+i-j];
			c2vi += dv2[j]*w2v[jp-jm][R+i-j];
			ct1 += dt1[j]*w3[jp-jm][R+i-j];
			ct2 += dt2[j]*(w2[jp-jm][R+i-j]-2.*w3[jp-jm][R+i-j]);
		}

		dep[i] = c3i + c2i - c2vi + ct1 + ct2 - mu;

		if (i>=R) dphi[i] = dep[i]+log(rho[i]/rho0);
		else dphi[i] = 0.0;

	}


}

double omega(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,om;
  double lR=R*dz;
  int i;
  
  weighted_densities();

  for (i=0;i<N-2*R;i++) {
    ln3=n3[i];
    ln2=n2[i];
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=n2v[i];
    lnv1=lnv2/(pi_4*lR);
    lnt1=nt1[i];
    lnt2=nt2[i];

    phi1=myphi1(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi2=myphi2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi3=myphi3(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,A,B);

    phi[i]=(phi1+phi2+phi3)+mode*p_bulk;

    if (i>=R) phiid[i]=rho[i]*(log(rho[i]/rho0)-1.0-mu);
    else phiid[i] = 0.0;

  }

  om = integrate_array_section(phi,dz,0,N-3*R)+integrate_array_section(phiid,dz,R,N-3*R);
  
	return om;
}

void minimize()
{
  double da1,hw=0.0,hw_old=1.0,hw_eps,hr;
  int i,j=1;

  hw_eps=1.e-21;

  for (int i=R;i<N;i++) rho[i]=rho0;
  for (int i=0;i<R;i++) rho[i]=0.0;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<ITS)) {
      
    hw_old=hw;

    weighted_densities();

    derivative();

    for (i=R;i<N-2*R;i++) {
      hr=rho[i];
      rho[i]=(1-a1)*rho[i]+a1*rho0*exp(-dep[i]);
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    j++;
  }
}

void gradient_descent()
{
  double hw=0.0,hw_old=1.0,hw_eps,hr,a1_max;
  int i,j=1;

  hw_eps=1.e-21;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<ITS)) {
      
    hw_old=hw;

    weighted_densities();

    derivative();

    for (i=R;i<N-2*R;i++) {
      hr=rho[i];

      a1_max = rho[i]/get_magnitude(dphi);

      if (a1 > a1_max) a1 = 0.1*a1_max;

      rho[i]-=a1*dphi[i];
      
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    j++;
  }
}

double adsorption() {
    
    int i;
    double ads[N];

    for (i=R; i<N-2*R; i++) ads[i] = rho[i]-rho0;

    return integrate_array_section(ads, dz, R, N-3*R)-rho0;
}


void adsorption_th() {

    double deta=1.e-6,eta_mid = eta;
    int count=0;

    dgdeta=0.0;

    while (count<=1) {

        eta = eta_mid + pow(-1, count+1)*deta;
				param_init();
        minimize();
        dgdeta += pow(-1.,count+1)*omega(1)/(2.*deta);
        count++;

    }

    eta = eta_mid;
}

void theorems_run() {

	FILE *fq;

  double step=1.e-2, hc;
	char adfln[250];
	sprintf(adfln, "ct_ads");
	fq = fopen(adfln, "w");

	for (double n = step; n<0.42; n+=step) {
        
	  eta = n;
		param_init();
		if (eta<0.35) a1 = 0.1;
		else if (eta<0.40) a1 = 0.01;
		else if (eta<0.45) a1 = 0.001;
	  else a1 = 0.0005;

	  minimize();
    adsorption_th();
	  Gamma = adsorption();

	  fprintf(fq, "%.4f	%.14f   %.14f	%.14f	%.14f   %.14f\n", eta, rho[R], p_bulk, Gamma, -dgdeta*dedmu, 4.*omega(1));

	}

	fclose(fq);

}

int main(int argc, char** argv) {

	void filename (char *fln);

  char fln[250],fname[50];
	double hc;

  if (argc >= 4) {
		eta = atof(argv[1]);
		A = atof(argv[2]);
		B = atof(argv[3]);
		if (argc >= 5) a1 = atof(argv[4]);
		else a1 = 0.1;

		filename(fname);

		init();
		minimize();
    
    Gamma = adsorption();
    adsorption_th();

		printf("Contact %.8f %.8f %.3e\n", rho[R], p_bulk, rho[R]-p_bulk);
    printf("Adsorption %.8f %.8f %.3e\n", Gamma, -dgdeta*dedmu, Gamma+dgdeta*dedmu);

		sprintf(fln, "hw_%s_%.3f_%.3f_%.3f",fname,eta,A,B);
		fp = fopen(fln, "w");

		for (int i=0; i<N-2*R; i++) {
		fprintf(fp, "%10.5lf %.14e %.14e %.14e %.14e\n", i*dz/2., 8.*rho[i], n3[i], dep[i], dphi[i]);
    }
    fclose(fp);

    } else {
			printf("\033[0;31m");
			printf("Please insert packing fraction and mixing parameter\n");
    }

    return 0;
}
