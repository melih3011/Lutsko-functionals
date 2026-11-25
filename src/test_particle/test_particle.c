#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8192
#define R 100
#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */

double rho[N],n3[N],n2[N],nv2[N],nt1[N],nt2[N];
double phi[N],phiid[N],dphi[N],dep[N],ads[N],Vext[N];
double w3[2*R+1][2*R+1],w2[2*R+1][2*R+1],wv2[2*R+1][2*R+1];
double wt21[2*R+1][2*R+1],wt22[2*R+1][2*R+1],wt23[2*R+1][2*R+1];
double eta,rho0,mu,xt,dr=1.0/(1.0*R),a1,p_bulk,Rb;
double A,B, rhoc, dgdeta, dedmu, Gamma;
int Rbi;

double get_magnitude(const double* arr) {

  double max = fabs(arr[0]);

  for (int i=1; i<N-2*R; ++i) {

    if (fabs(arr[i]) > max) max = fabs(arr[i]);
  }

  return max;

}

double fr2int(const double *f,double ldr,int n, double a)
{
  double val=0.0,hr1=6.*a*a,hr2=4*a*ldr,hr3=ldr*ldr;

  for (int j=0;j<n;j++) {
    val += f[j]*(hr1+hr2*(1+3*j)+hr3*(1+4*j+6*j*j))+ 
    f[j+1]*(hr1+hr2*(2+3*j)+hr3*(3+8*j+6*j*j));
  }
  
  return ldr*val/12.0;
}


void init_matrix(const int n, double wm[2*n+1][2*n+1], const double *lwm) {
  
  double p1[2], p2[3], p3[4], p4[5];
  p1[0] = 0.5; p1[1] = 0.5;
  p2[0] = 1./3.; p2[1] = 4./3.; p2[2] = 1./3.;
  p3[0] = 3./8.; p3[1] = 9./8; p3[2] = 9./8.; p3[3] = 3./8.;
  p4[0] = 14./45.; p4[1] = 64./45.; p4[2] = 24./45.; p4[3] = 64./45.; p4[4] = 14./45.;
  
  
  for (int i=0;i<5;++i) {
    for (int j=0;j<2*n;++j) {
      wm[i][j]=0.0;
    }
  }
  
  for (int j=0; j<2; j++) wm[1][j] = lwm[j]*p1[j];
  for (int j=0; j<3; j++) wm[2][j] = lwm[j]*p2[j];
  for (int j=0; j<4; j++) wm[3][j] = lwm[j]*p3[j];
  for (int j=0; j<5; j++) wm[4][j] = lwm[j]*p4[j];
  
  for (int i=5; i<=2*n; i++) {
    for (int j=0; j<2*n-i; j++) {
      wm[i][2*n-j] = 0.0;
    }
    
    wm[i][i] = 3./8.*lwm[i]; wm[i][i-1] = 7./6.*lwm[i-1]; wm[i][i-2] = 23./24.*lwm[i-2]; 
    
    for (int j=3; j<=i-3; j++) wm[i][j] = lwm[j];
    
    wm[i][2] = 23./24.*lwm[2]; wm[i][1] = 7./6.*lwm[1]; wm[i][0] = 3./8.*lwm[0];
  }
  
}

double sk0()
{
  double hr[N];

  for (int i=R;i<=N-2*R;i++) hr[i-R]=rho[i]-rho0;

  double sum=fr2int(hr,dr,N-4*R,Rb+1.)-pow(Rb+1.,3)*rho0/3.;

  return 1.0+pi_4*sum;
}

void param_init() {

  double pressure(double eta,double rho0,double a,double b);
  double chemical_potential(double eta,double a,double b);
  double sk0_bulk(double eta,double a,double b);
  double detadmu(double eta, double a, double b);
  
  rho0=3.0*eta/pi_4;
  mu=chemical_potential(eta,A,B);
  xt=sk0_bulk(eta,A,B);
  p_bulk=pressure(eta,rho0,A,B);
  dedmu=detadmu(eta,A,B);
}


void data_init()
{
  double lw3[2*R+1],lw2[2*R+1],lwv2[2*R+1];
  double lwt21[2*R+1],lwt22[2*R+1],lwt23[2*R+1];
  double lR=R*dr;

  param_init();

  for (int i=0;i<N;i++) {
    n3[i]=0.0;
    n2[i]=0.0;
    nv2[i]=0.0;
    nt1[i]=0.0;
    nt2[i]=0.0;
    Vext[i]=0.0;
  }
  
  for (int i=0;i<=2*R;i++) {
    double z=(i-R)*dr;
    lw3[i]=dr*pi*(lR*lR-z*z);
    lw2[i]=dr*pi_2*lR;
    lwv2[i]=dr*pi_2*z;
    //tensor contributions
    lwt21[i]=dr*pi_2*0.25*pow(z*z-lR*lR,2)/lR;
    lwt22[i]=-dr*pi_2*z*(z*z-lR*lR)/lR;
    lwt23[i]=dr*pi_2*z*z/lR;
  }

  init_matrix(R, w3, lw3);
  init_matrix(R, w2, lw2);
  init_matrix(R, wv2, lwv2);
  init_matrix(R, wt21, lwt21);
  init_matrix(R, wt22, lwt22);
  init_matrix(R, wt23, lwt23);

}

void weighted_densities()
{
  n3[0]= 0.0;
  n2[0]= 0.0;
  nv2[0]= 0.0;
  nt1[0]=0.0;
  nt2[0]=0.0;
  
  for (int i = 1; i < N-R ; i++) {

    double hr = Rb+i*dr, hrho;
    
    double n3i = 0.0, n2i = 0.0;
    double nv2i = 0.0, nt2i = 0.0;

    int jp=i+R;
    int jm=i-R;
    if (jm<R) jm=R;
    
    for (int j = jm; j <= jp; j++) {
      hrho = rho[j]*(Rb+j*dr);
      n3i += hrho*w3[jp-jm][R+i-j];
      n2i += hrho*w2[jp-jm][R+i-j];
      nv2i += hrho*wv2[jp-jm][R+i-j];
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

  double d2[N],d3[N],dv2[N],dt1[N],dt2[N];
  double lR=R*dr;
  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;

  for (int i=0;i<N-R;i++) {
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
  for (int i=N-R;i<N;i++) {
    d3[i]=0.0;
    d2[i]=0.0;
    dv2[i]=0.0;
    dt1[i]=0.0;
    dt2[i]=0.0;
  }

  for (int i = 0; i < N-R; i++) {

    double hr = Rb+i*dr;
    double c3i = 0.0, c2i = 0.0;
    double cv21i = 0.0, cv22i = 0.0;
    double ct1 = 0.0, ct2 = 0.0;
    
    int jp = i + R;
    int jm = i - R;
    if (jm<0) jm=0;
    
    for (int j = jm; j <= jp; j++) {
      double hrs = Rb+j*dr;

      c3i += d3[j]*hrs*w3[jp-jm][(R+i)-j];
      c2i += d2[j]*hrs*w2[jp-jm][(R+i)-j];
      cv21i += dv2[j]*w3[jp-jm][(R+i)-j]; 
      cv22i += dv2[j]*hrs*wv2[jp-jm][(R+i)-j];
      ct1 += dt1[j]*hrs*0.5*w2[jp-jm][R+i-j]-dt1[j]*0.5*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
      ct2 += dt2[j]*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
    }
    
    dep[i] = (c3i + c2i - cv22i + cv21i + ct1 + ct2)/hr - mu;
    
    if (i>=R) dphi[i]=dep[i]+log(rho[i]/rho0);
    else dphi[i]=0.0;
  } 
}

double surface_tension(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,om;
  double lR=R*dr;

  for (int i=0;i<N-2*R;i++) {
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

  om = fr2int(phi,dr,N-3*R-1,Rb)+fr2int(phiid,dr,N-4*R-1,Rb+1.);
  
  if (mode==0) return pi_4*om;

  return om/pow(Rb+1.,2);
}

double omega(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,om;
  double lR=R*dr;

  for (int i=0;i<N-2*R;i++) {
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

  om = fr2int(phi,dr,N-3*R-1,Rb)+fr2int(phiid,dr,N-4*R-1,Rb+1.);
  
  if (mode==0) return pi_4*om;

  return om/pow(Rb,2);
}

double adsorption(int mode)
{

  for (int i=R;i<=N-2*R;i++) ads[i-R]=rho[i]-mode*rho0;

  return fr2int(ads,dr,N-4*R-1,Rb+1.);
}


void gradient_descent()
{
  double hw=0.0,hw_old=1.0,hw_eps,hr;
  int j=1;

  hw_eps=1.e-21;

  for (int i=R;i<N;i++) rho[i]=rho0;
  for (int i=0;i<R;i++) rho[i]=0.0;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<50000)) {

    hw_old=hw;
    weighted_densities();
    derivative();

    for (int i=R;i<N-2*R;i++) {
      hr=rho[i];
      double a1_max=rho[i]/get_magnitude(dphi);
      if (a1 > a1_max) a1 = 0.5*a1_max;
      rho[i]-=a1*dphi[i];
      
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    j++;
  }
}


void minimize()
{
  double hw=0.0,hw_old=1.0,hw_eps=1.e-18,hr;
  int j = 1, max_its=20000;

  for (int i=R;i<N;i++) rho[i]=rho0;
  for (int i=0;i<R;i++) rho[i]=0.0;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<max_its)) {

    hw_old=hw;

    weighted_densities();

    derivative();

    for (int i=R;i<N-2*R;i++) {
      hr=rho[i];
      rho[i]=(1-a1)*rho[i]+a1*rho0*exp(-dep[i]);
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    j++;

  }
  
  if (j>=max_its) printf("Exceeded maximum number of iterations!\n");
}


void adsorption_th() {

  double deta=1.e-6, etamid=eta;
  int count=0;
  
  dgdeta = 0.0;

  while (count<=1) {
    eta = etamid -pow(-1,count)*deta;
    param_init();
    minimize();
    dgdeta += pow(-1.,count+1)*surface_tension(1)/(2.*deta);
    count++;
  }

  eta = etamid;
}

void contact_th() {

  double dgdR=0.0, lR=R*dr, Rmid=Rb, dR=1.e-6, hgamma=surface_tension(1);
  int count=0;
  
  while (count<=1) {
    Rb = Rmid - pow(-1,count)*dR;
    param_init();
    minimize();
    dgdR += pow(-1.,count+1)*surface_tension(1)/(2.*dR);
    count++;
  }

  Rb = Rmid;
  rhoc = p_bulk+2.*hgamma/(Rb+lR)+dgdR;
}

int main(int argc, char **argv)
{
  void filename(char *fln);
  char fname[50];

  filename(fname);

  if (argc>=6) {

    Rb=atof(argv[1]);
    eta=atof(argv[2]);
    A=atof(argv[3]);
    B=atof(argv[4]);
    int choice=atoi(argv[5]);
    
    if (argc>=7) a1=atof(argv[6]);
    else a1=0.1;

    Rbi=(int)(R*Rb);

    switch (choice) {
      default: {

        FILE *fp;
        char fln[250];

        data_init();
        minimize();
        //gradient_descent();

        sprintf(fln,"sp_%s_%6.4lf_%6.4lf",fname,eta,Rb);
        fp=fopen(fln,"w");
        fprintf(fp,"# Rb       = %10.8lf\n",Rb);
        fprintf(fp,"# p_bulk   = %10.8lf\n",p_bulk);
        fprintf(fp,"# mu_ex,bl = %10.8lf\n",mu);
        fprintf(fp,"# mu_ex,om = %10.8lf\n",omega(1)*pi_4 + pi_4*p_bulk/3.);
        fprintf(fp,"# Sk(0)_bl = %10.8lf\n",xt);
        fprintf(fp,"# Sk(0)_om = %10.8lf\n",sk0());
        fprintf(fp,"# Gamma = %10.8lf\n",adsorption(1));

        for (int i=0;i<N-2*R;i++) {
          fprintf(fp,"%.14e %.14e %.14e %.14e %.14e\n",(Rbi+i)*dr/2.,rho[i]/rho0,n3[i],phi[i],dep[i]);
        }
        fclose(fp);

        Gamma = adsorption(1)/pow(Rb+1.,2);
        contact_th();
        adsorption_th();

        printf("Accuracy %.2e\n", get_magnitude(dphi));
        printf("Contact %.6f %.6f %.3e\n", rho[R], rhoc, fabs(rho[R]/rhoc-1.));
        printf("Adsorption %.6f %.6f %.3e\n", Gamma, -dgdeta*dedmu, fabs(-Gamma/(dgdeta*dedmu)-1.));
        printf("Excess chemical potential %.6f %.6f\n", mu, omega(1)*pi_4 + pi_4*p_bulk/3.);
        printf("Isothermal Compressibility %.6f %.6f\n", xt, sk0());

        break;
      }

      case 1: {

        FILE *frf;
        char rf[250];

        data_init();

        sprintf(rf,"%s_run_eta",fname);
        frf=fopen(rf,"w");

        for (int n=1; n<(int)(0.49/0.01); n++) {

          eta = n*0.01;
          param_init();

          if (eta < 0.35) a1 = 0.1;
          else if (eta < 0.40) a1 = 0.01;
          else if (eta < 0.45) a1 = 0.005;
          else a1 = 0.001;

          minimize();

          double mu_dev = (omega(1)*pi_4+pi_4*p_bulk/3.)/mu-1.;
          double xt_dev = sk0()/xt-1.;

          printf("%.4f %.14e %.14e, Accuracy %.1e\n", eta, mu_dev, xt_dev, get_magnitude(dphi));
          fprintf(frf,"%.4f %.14e %.14e %.14e %.14e %.14e %.14e\n", eta, omega(1)*pi_4 + pi_4*p_bulk/3., mu, mu_dev, sk0(), xt, xt_dev);
        }
        fclose(frf);
        break;
      }

      case 2: {

        double Am=1.0, Ap=1.4, dA=5.e-2, Bm=-1.5, Bp=1.5;
        int stepsA=(int)((Ap-Am)/dA), stepsB=(int)((Bp-Bm)/dA);

        data_init();

        for (int k=0; k<=stepsA; ++k) {
          A = Am+k*dA;
          //Bm = (9.-8.*A)/2.;
          //for (double b=Bm-0.3; b<=Bm+0.3; b+=dA) {
          for (int l=0; l<=stepsB; ++l) {

            B = Bm+l*dA;
            //B = b;
            param_init();
            minimize();
            printf("%.3f  %.3f    %.14e  %.14e\n", A, B, ((omega(1)*pi_4 + pi_4*p_bulk/3.)/mu-1.), (sk0()/xt-1.));
          }
          printf(" \n");
        }
        break;
      }
    }

  } else fprintf(stderr,"syntax: %s Rb eta [a]\n",argv[0]);
}


