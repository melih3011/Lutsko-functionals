#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 10000//4096
#define R 50

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */

double rho[N],n3[N],n2[N],nv2[N],nt1[N],nt2[N];
double phi[N],phiid[N],dphi[N],dep[N],ads[N],Vext[N];
double w3[2*R+1][2*R+1],w2[2*R+1][2*R+1],wv2[2*R+1][2*R+1];
//tensorial contributions must be separated into several terms exhibiting different 1/r-dependencies
//it holds that t1 = 1/2 (n2-t2), hence former is fixed by latter weighted densities
double wt21[2*R+1][2*R+1],wt22[2*R+1][2*R+1],wt23[2*R+1][2*R+1];
double eta,rho0,mu,dr=1.0/(1.0*R),a1,p_bulk,Rb;
int Rbi;

double A,B;
double rhoc;
double dgdeta, dedmu, Gamma;

FILE *fp, *frf, *frun;

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

double integrate(const double *f, double dx, int i0, int i1) {

	double sum=0.0;
	int i;

	if (i0 == i1) return 0.0;
	else {
	
		//Eq.4.1.18 of Numerical Recipes

		sum=55./24.*(f[i0+1]+f[i1-1])-1./6.*(f[i0+2]+f[i1-2])+11./8.*(f[i0+3]+f[i1-3]);

		for (i=i0+4; i<=i1-4; ++i) sum+=f[i];

		return sum*dx;

	}

}

double fr2int(double f[],double ldr,int n, double a)
{
  double val=0.0,hr1=6*a*a,hr2=4*a*ldr,hr3=ldr*ldr;
  register int j;

  for (j=0;j<n;j++)
    val+=f[j]*(hr1+hr2*(1+3*j)+hr3*(1+4*j+6*j*j))+
      f[j+1]*(hr1+hr2*(2+3*j)+hr3*(3+8*j+6*j*j));
  return(ldr*val/12.0);
}


double fr2int_2(const double *f,double ldr,int n,double a) {

    double val=0.0, w[4], hr1, hr2;
    int i;

    w[0]=17./48., w[1]=59./48., w[2]=43./48., w[3]=49./48.;

    
    for (i=0;i<=n-4;++i) {
        
        hr1 = pow(a+i*ldr,2);
        
        if (i<4) {
            hr2 = pow(a+(n-i)*ldr,2);
            val+=w[i]*(f[i]*hr1+f[n-i]*hr2);
        } else val+=f[i]*hr1;
    }
	
	return val*ldr;

}

double sk0()
{
  double hr[N],sum;
  int i;

  for (i=R;i<N-3*R;i++) hr[i-R]=rho[i]-rho0;
    
  sum=fr2int_2(hr,dr,N-4*R-1,Rb+1.)-pow(Rb+1.,3)*rho0/3.;

  return 1.0+pi_4*sum;
}

void data_init()
{

  double pressure(double eta,double rho0,double a,double b);
  double chemical_potential(double eta,double a,double b);
  double detadmu(double eta, double a, double b);

  double lw3[2*R+1],lw2[2*R+1],lwv2[2*R+1];
  double lwt21[2*R+1],lwt22[2*R+1],lwt23[2*R+1];
  int i,j;
  double z,lR=R*dr;

  rho0=3.0*eta/(pi_4);

  mu=chemical_potential(eta,A,B);

  p_bulk=pressure(eta,rho0,A,B);

  for (i=0;i<N;i++) {
    n3[i]=0.0;
    n2[i]=0.0;
    nv2[i]=0.0;
	nt1[i]=0.0;
	nt2[i]=0.0;
    Vext[i]=0.0;
  }
  
  for (i=0;i<=2*R;i++) {
    z=(i-R)*dr;
    lw3[i]=dr*pi*(lR*lR-z*z);
    lw2[i]=dr*pi_2*lR;
    lwv2[i]=dr*pi_2*z;
    //tensor contributions
    lwt21[i]=dr*pi_2*0.25*pow(z*z-lR*lR,2)/lR;
    lwt22[i]=-dr*pi_2*z*(z*z-lR*lR)/lR;
    lwt23[i]=dr*pi_2*z*z/lR;
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
			w3[i][j] = 0.0; w2[i][j] = 0.0; wv2[i][j] = 0.0;
			wt21[i][j] = 0.0; wt22[i][j] = 0.0; wt23[i][j] = 0.0;
		}
	}
	

  for (int j=0; j<2; j++) {
    w3[1][j] = lw3[1]*p1[j]; w2[1][j] = lw2[1]*p1[j]; wv2[1][j] = lwv2[1]*p1[j];
    wt21[1][j] = lwt21[1]*p1[j]; wt22[1][j] = lwt22[1]*p1[j]; wt23[1][j] = lwt23[1]*p1[j];
  }

  for (int j=0; j<3; j++) {
    w3[2][j] = lw3[2]*p2[j]; w2[2][j] = lw2[2]*p2[j]; wv2[2][j] = lwv2[2]*p2[j];
    wt21[2][j] = lwt21[2]*p2[j]; wt22[2][j] = lwt22[2]*p2[j]; wt23[2][j] = lwt23[2]*p2[j];
  }

  for (int j=0; j<4; j++) {
    w3[3][j] = lw3[3]*p3[j]; w2[3][j] = lw2[3]*p3[j]; wv2[3][j] = lwv2[3]*p3[j];
    wt21[3][j] = lwt21[3]*p3[j]; wt22[3][j] = lwt22[3]*p3[j]; wt23[3][j] = lwt23[3]*p3[j];
  }

  for (int j=0; j<5; j++) {
    w3[4][j] = lw3[4]*p4[j]; w2[4][j] = lw2[4]*p4[j]; wv2[4][j] = lwv2[4]*p4[j];
    wt21[4][j] = lwt21[4]*p4[j]; wt22[4][j] = lwt22[4]*p4[j]; wt23[4][j] = lwt23[4]*p4[j];
  }

	for (int i=5; i<=2*R; i++) {
		for (int j=0; j<2*R-i; j++) {
			w3[i][2*R-j] = 0.0; w2[i][2*R-j] = 0.0; wv2[i][2*R-j] = 0.0;
			wt21[i][2*R-j] = 0.0; wt22[i][2*R-j] = 0.0; wt23[i][2*R-j] = 0.0;
		}
	
	w3[i][i] = 3./8.*lw3[i]; w3[i][i-1] = 7./6.*lw3[i-1]; w3[i][i-2] = 23./24.*lw3[i-2]; 
	w2[i][i] = 3./8.*lw2[i]; w2[i][i-1] = 7./6.*lw2[i-1]; w2[i][i-2] = 23./24.*lw2[i-2]; 
	wv2[i][i] = 3./8.*lwv2[i]; wv2[i][i-1] = 7./6.*lwv2[i-1]; wv2[i][i-2] = 23./24.*lwv2[i-2]; 

	wt21[i][i] = 3./8.*lwt21[i]; wt21[i][i-1] = 7./6.*lwt21[i-1]; wt21[i][i-2] = 23./24.*lwt21[i-2];
	wt22[i][i] = 3./8.*lwt22[i]; wt22[i][i-1] = 7./6.*lwt22[i-1]; wt22[i][i-2] = 23./24.*lwt22[i-2];
	wt23[i][i] = 3./8.*lwt23[i]; wt23[i][i-1] = 7./6.*lwt23[i-1]; wt23[i][i-2] = 23./24.*lwt23[i-2];

	for (int j=3; j<=i-3; j++) {
		w3[i][j] = lw3[j]; w2[i][j] = lw2[j]; wv2[i][j] = lwv2[j];
		wt21[i][j] = lwt21[j]; wt22[i][j] = lwt22[j]; wt23[i][j] = lwt23[j];
	} 

	w3[i][2] = 23./24.*lw3[2]; w3[i][1] = 7./6.*lw3[1]; w3[i][0] = 3./8.*lw3[0];
	w2[i][2] = 23./24.*lw2[2]; w2[i][1] = 7./6.*lw2[1]; w2[i][0] = 3./8.*lw2[0];
	wv2[i][2] = 23./24.*lwv2[2]; wv2[i][1] = 7./6.*lwv2[1]; wv2[i][0] = 3./8.*lwv2[0];
	
  wt21[i][2] = 23./24.*lwt21[2]; wt21[i][1] = 7./6.*lwt21[1]; wt21[i][0] = 3./8.*lwt21[0];
  wt22[i][2] = 23./24.*lwt22[2]; wt22[i][1] = 7./6.*lwt22[1]; wt22[i][0] = 3./8.*lwt22[0];
  wt23[i][2] = 23./24.*lwt23[2]; wt23[i][1] = 7./6.*lwt23[1]; wt23[i][0] = 3./8.*lwt23[0];

	}

}

void weighted_densities()
{
  int i, j, jp, jm;
  double hr,hrho,n3i,n2i,nv2i,nt1i,nt2i;

  n3[0]= 0.0;
  n2[0]= 0.0;
  nv2[0]= 0.0;
  nt1[0]=0.0;
  nt2[0]=0.0;
  
  for(i = 1; i < N-R ; i++) {

    hr = (Rbi+i)*dr;
    
    n3i = 0.;
    n2i = 0.;
    nv2i = 0.;
    nt1i = 0.;
    nt2i = 0.;

    jp=i+R;
    jm=i-R;
    if (jm<R) jm=R;
    
    for(j = jm; j <= jp; j++) {
      hrho = rho[j]*(Rbi+j)*dr;
      n3i += hrho*w3[jp-jm][R+i-j];
      n2i += hrho*w2[jp-jm][R+i-j];
      nv2i += hrho*wv2[jp-jm][R+i-j];
      //tensor contribution where only t2 has to be computed
      nt2i += hrho*(wt21[jp-jm][R+i-j]/pow(hr,3)+wt22[jp-jm][R+i-j]/pow(hr,2)+wt23[jp-jm][R+i-j]/hr);
    }
    
    n3[i] = n3i/hr;
    n2[i] = n2i/hr;
    nv2[i] = (n3[i]+nv2i)/hr;
    nt1[i] = 0.5*(n2i/hr-nt2i);
    nt2[i] = nt2i;
  }
}

void derivative()
{
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
  double lR=R*dr,hr,hrs;
  double c3i,c2i,cv21i,cv22i,ct1,ct2;
  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;

  for (i=0;i<N-R;i++) {
    ln3=n3[i];
    ln2=n2[i];      
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);
    lnt1=nt1[i];
    lnt2=nt2[i];  

    d2[i]=myd2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,lR,A,B);
    d3[i]=myd3(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,lR,A,B);
    dv2[i]=mydv2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,lR,A,B);
    dt1[i]=mydt1(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,lR,A,B);
    dt2[i]=mydt2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,lR,A,B);
  }
  for (i=N-R;i<N;i++) {
    d3[i]=0.0;
    d2[i]=0.0;
    dv2[i]=0.0;
    dt1[i]=0.0;
    dt2[i]=0.0;
  }

  for(i = 0; i < N-R; i++) {  

    hr = (Rbi+i)*dr;
    c3i = 0.;
    c2i = 0.;
    cv21i = 0.;
    cv22i = 0.;
    ct1 = 0.0;
    ct2 = 0.0;
    
    jp = i + R;
    jm = i - R;
    if (jm<0) jm=0;
    
    for(j = jm; j <= jp; j++) {
      hrs = (Rbi+j)*dr;

      c3i += d3[j]*hrs*w3[jp-jm][(R+i)-j];
      c2i += d2[j]*hrs*w2[jp-jm][(R+i)-j];
      cv21i += dv2[j]*w3[jp-jm][(R+i)-j]; //i <-> j
      cv22i += dv2[j]*hrs*wv2[jp-jm][(R+i)-j]; //i <-> j
      ct1 += dt1[j]*hrs*0.5*w2[jp-jm][R+i-j]-dt1[j]*0.5*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
      ct2 += dt2[j]*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
    }
    
    dep[i] = (c3i + c2i - cv22i + cv21i + ct1 + ct2)/hr - mu;
    
    if (i>=R) dphi[i]=dep[i]+log(rho[i]/rho0);
    else dphi[i]=0.0;
  } 
}

//additional method to compute surface tension gamma with right definitons at the boundary
double gamma(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double in3,in2,inv2,hr,ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,om;
  double lR=R*dr;
  int i;
  
  weighted_densities();

  for (i=0;i<N-2*R;i++) {
    ln3=n3[i];
    ln2=n2[i];
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);
    lnt1=nt1[i];
    lnt2=nt2[i];

    phi1=myphi1(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi2=myphi2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi3=myphi3(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,A,B);

    phi[i]=(phi1+phi2+phi3);

    if (i>=R) phiid[i-R]=rho[i]*(log(rho[i]/rho0)-1.0-mu)+mode*p_bulk;

  }

  om = fr2int_2(phi,dr,N-3*R-1,Rb)+fr2int_2(phiid,dr,N-4*R-1,Rb+1.);
  
  if (mode==0) return pi_4*om;
  else return om/pow(Rb+lR,2);
}

double omega(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double in3,in2,inv2,hr,ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,om;
  double lR=R*dr;
  int i;
  
  weighted_densities();

  for (i=0;i<N-2*R;i++) {
    ln3=n3[i];
    ln2=n2[i];
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);
    lnt1=nt1[i];
    lnt2=nt2[i];

    phi1=myphi1(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi2=myphi2(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2);
    phi3=myphi3(ln0,ln1,ln2,ln3,lnv1,lnv2,lnt1,lnt2,A,B);

    phi[i]=(phi1+phi2+phi3)+mode*p_bulk;

    if (i>=R) phiid[i-R]=rho[i]*(log(rho[i]/rho0)-1.0-mu);

  }

  om = fr2int_2(phi,dr,N-3*R-1,Rb)+fr2int_2(phiid,dr,N-4*R-1,Rb+1.);
  
  if (mode==0) return pi_4*om;
  else return om/pow(Rb,2);
}

double adsorption(int mode)
{
  int i;

  for (i=R;i<=N-2*R;i++) ads[i-R]=rho[i]-mode*rho0;

  return fr2int_2(ads,dr,N-3*R,Rb+1.);
}

void minimize()
{
  double da1,hw=0.0,hw_old=1.0,hw_eps,hr;
  int i,j=1;

  hw_eps=0.000000000000000001;

  for (i=R;i<N;i++) rho[i]=rho0;
  for (i=0;i<R;i++) rho[i]=0.0;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<10000)) {
  //  for (j=0;j<500;j++) {
    hw_old=hw;

    // hw=omega(0);

    weighted_densities();

    derivative();

    for (i=R;i<N-4*R;i++) {
      hr=rho[i];
      rho[i]=(1-a1)*rho[i]+a1*rho0*exp(-dep[i]);
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    j++;
  }
}


void contact_th() {

  double dgdR=0., lR=R*dr, dR=dr;
  double hgamma=gamma(1);

  int count=0;
  while (count<=1) {

    Rb += -dR + 3.*count*dR;
    Rbi = (int)(Rb*R);

    minimize();

    dgdR += pow(-1.,count+1)*gamma(1)/(2.*dR);
    count++;
  }

  Rb -= dR;
  Rbi = (int)(Rb*R);

  rhoc = p_bulk+2.*hgamma/(Rb+lR)+dgdR;

}


void adsorption_th() {

  double deta=1.e-6;
	double old_eta=eta;
	int count=0;
  
    dedmu = detadmu(eta, A, B);
    dgdeta = 0.;

	while (count<=1) {

        eta += -deta + 3.*count*deta;

		rho0 = 3.*eta/pi_4;
		mu = chemical_potential(eta, A, B);
		p_bulk = pressure(eta, rho0, A, B);

		minimize();

		dgdeta += pow(-1.,count+1)*gamma(1)/(2.*deta);

		count++;
	}

	eta = old_eta;

}


void theorems_run () {

  FILE *fq;
  char adfln[250];

  double step=1.e-2;

  sprintf(adfln, "ct_ads");
  fq = fopen(adfln, "w");

  for (double n = 0.3; n<0.4; n+=step) {

    if (n < 0.35) a1 = 0.1;

    else if (n < 0.40) a1 = 0.05;

    else a1 = 0.001;

    eta = n;
    rho0 = 3.*eta/pi_4;
    mu = chemical_potential(eta, A, B);
    p_bulk = pressure(eta, rho0, A, B);

    dedmu = detadmu(eta, A, B);

    minimize();

    Gamma = adsorption(1)/pow(Rb+1.,2);

    contact_th();
    adsorption_th();
		
    printf("%.4f    %.14e    %.14e	%.14e	%.14e\n", eta, rho[R], rhoc, Gamma, -dgdeta*dedmu);
    fprintf(fq, "%.4e	%.14e    %.14e    %.14e	%.14e\n", eta, rho[R], rhoc, Gamma, -dgdeta*dedmu);
    //fprintf(fq, "%.4e	%.14e    %.14e    %.14e	%.14e\n", eta, p_bulk*8., rho0*8.*(1.+eta+eta*eta)/pow(1.-eta,3), rho0*8.*(1.+eta+eta*eta-eta*eta*eta)/pow(1.-eta,3));
  
  }

  fclose(fq);

}


int main(int argc,char *argv[])
{
  int choice;
  void filename(char *fln);
  double sk0_bulk(double eta,double A,double B);
  
  
  int i,j;
  char fln[250],fname[50],rf[250],run[250],command[250];
  double r,hw;

  double Am=-2., Ap=2., dA=1.e-2;
  int steps = (int)((Ap-Am)/dA), counter=0;

  if (argc>=6) {
    Rb=atof(argv[1]);
    eta=atof(argv[2]);
    A = atof(argv[3]), B = atof(argv[4]);
    choice = atoi(argv[5]);
    if (argc>=7) a1=atof(argv[6]);
    else a1=0.1;

    Rbi=(int)(R*Rb);

    data_init();
    filename(fname);

    switch (choice) {

    case 0:

      minimize();
      

      sprintf(fln,"sp_%s_%6.4lf_%6.4lf",fname,eta,Rb);
      fp=fopen(fln,"w");
      fprintf(fp,"# Rb       = %10.8lf\n",Rb);
      fprintf(fp,"# p_bulk   = %10.8lf\n",p_bulk);
      fprintf(fp,"# mu_ex,bl = %10.8lf\n",mu);
      fprintf(fp,"# mu_ex,om = %10.8lf\n",omega(1)*pi_4 + pi_4*p_bulk/3.);
      fprintf(fp,"# Sk(0)_bl = %10.8lf\n",sk0_bulk(eta,A,B));
      fprintf(fp,"# Sk(0)_om = %10.8lf\n",sk0());
      fprintf(fp,"# Gamma = %10.8lf\n",adsorption(1));

      for (i=0;i<N-5*R;i++) {
        fprintf(fp,"%15.10lf   %20.19e  %20.19e  %20.19lg %20.19lg %10.19lg %lg %lg  %lg\n",(Rbi+i)*dr/2.,rho[i]*8.,rho[i]/rho0,n3[i],n2[i],nv2[i],nt2[i],dphi[i], dep[i]);
      }
      fclose(fp);
        
      sprintf(command,"mv %s ./data",fln);  
      system(command);

  		printf("%.14e	%.14e\n", fabs((omega(1)*pi_4 + pi_4*p_bulk/3.)/mu-1.)*100., fabs(sk0()/sk0_bulk(eta,A,B)-1.)*100.);
      break;

    case 1:

      theorems_run();
      break;

    case 2: 

      sprintf(rf,"lk_run_eta");
      frf=fopen(rf,"w");

      for(double n=0.; n<0.45; n+=1.e-3) {
          eta = n;
          rho0=3.0*eta/(pi_4);
          mu=chemical_potential(eta,A,B);
          p_bulk=pressure(eta,rho0,A,B);

          if (eta > 0.4) a1 = 0.05;

          minimize();

          fprintf(frf,"%.4f   %.14e   %.14e %.14e    %.14e\n", eta, omega(1)*pi_4 + pi_4*p_bulk/3., mu, sk0(), sk0_bulk(eta,A,B));
      }
      fclose(frf);
      break;


    case 3:

      printf("Calculating at %.2f\n", eta);

      sprintf(run,"run_%s_%6.4lf_%6.4lf",fname,eta,Rb);
      frun=fopen(run,"w");

      fprintf(frun, "#A dmu dchit\n");
      for (int k=0; k<=steps; ++k) {
      
        A = Am+k*dA;
      
        for (int l=0; l<=steps; ++l) {

          B = Am+l*dA;
          //bulk pressure
          p_bulk = pressure(eta,rho0,A,B);
          //chemcial potential
          mu = chemical_potential(eta,A,B);

          minimize();

          fprintf(frun,"%.3f  %.3f    %.14e  %.14e\n", A, B, ((omega(1)*pi_4 + pi_4*p_bulk/3.)/mu-1.)*100., (sk0()/sk0_bulk(eta,A,B)-1.)*100.);

          if (l == steps) fprintf(frun, " \n");

          counter++;
          printf("\rProgress %.3f pc", counter/pow(steps+1,2)*100.);
          fflush(stdout);

        }

      
      }

      fclose(frun);
      printf("\n");
      break;
	

  }	

  } else fprintf(stderr,"syntax: %s Rb eta [a]\n",argv[0]);
}


