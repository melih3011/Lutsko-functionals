#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */

#define R 100

double *rho,*n3,*n2,*nv2,*nt1,*nt2;
double *phi,*phiid,*dphi,*dep,*ads,*Vext;
double w3[2*R+1][2*R+1],w2[2*R+1][2*R+1],wv2[2*R+1][2*R+1];
double wt21[2*R+1][2*R+1],wt22[2*R+1][2*R+1],wt23[2*R+1][2*R+1];
double eta,rho0,mu,xt,dr=1.0/(1.0*R),a1,p_bulk;
double Nr,Rb,A,B;
int N,Rbi;

FILE *fp;

double fr2int(const double *f, double ldr, int n, double a)
{
  double val=0.0,hr1=6*a*a,hr2=4*a*ldr,hr3=ldr*ldr;
  register int j;

  for (j=0;j<n;j++)
    val+=f[j]*(hr1+hr2*(1+3*j)+hr3*(1+4*j+6*j*j))+
      f[j+1]*(hr1+hr2*(2+3*j)+hr3*(3+8*j+6*j*j));
  return(ldr*val/12.0);
}


void param_init() {
 
  double pressure(double eta,double rho0,double a,double b);
  double chemical_potential(double eta,double a,double b);
  
  rho0=3.0*eta/(pi_4);
  mu=chemical_potential(eta,A,B);
  p_bulk=pressure(eta,rho0,A,B);
}

void init_matrix(int n, double wm[2*n+1][2*n+1], const double *lwm) {
  
  int i,j;
  double p1[2], p2[3], p3[4], p4[5];
  p1[0] = 0.5; p1[1] = 0.5;
  p2[0] = 1./3.; p2[1] = 4./3.; p2[2] = 1./3.;
  p3[0] = 3./8.; p3[1] = 9./8; p3[2] = 9./8.; p3[3] = 3./8.;
  p4[0] = 14./45.; p4[1] = 64./45.; p4[2] = 24./45.; p4[3] = 64./45.; p4[4] = 14./45.;
  
  
  for (i=0;i<5;++i) {
    for (j=0;j<2*n;++j) {
      wm[i][j]=0.0;
    }
  }
  
  for (j=0; j<2; j++) wm[1][j] = lwm[j]*p1[j];
  for (j=0; j<3; j++) wm[2][j] = lwm[j]*p2[j];
  for (j=0; j<4; j++) wm[3][j] = lwm[j]*p3[j];
  for (j=0; j<5; j++) wm[4][j] = lwm[j]*p4[j];
  
  for (i=5; i<=2*n; i++) {
    for (j=0; j<2*n-i; j++) {
      wm[i][2*n-j] = 0.0;
    }
    
    wm[i][i] = 3./8.*lwm[i]; wm[i][i-1] = 7./6.*lwm[i-1]; wm[i][i-2] = 23./24.*lwm[i-2]; 
    
    for (j=3; j<=i-3; j++) wm[i][j] = lwm[j];
    
    wm[i][2] = 23./24.*lwm[2]; wm[i][1] = 7./6.*lwm[1]; wm[i][0] = 3./8.*lwm[0];
  }
  
}

void data_init()
{
  double lw3[2*R+1],lw2[2*R+1],lwv2[2*R+1];
  double lwt21[2*R+1],lwt22[2*R+1],lwt23[2*R+1];
  int i,j;
  double z,lR=R*dr;
  
  N = (int)(Nr/dr);
  Rbi = N-R;
  Rb = Nr-1.;

  rho = (double*)malloc(N*sizeof(double));
  n3 = (double*)malloc(N*sizeof(double));
  n2 = (double*)malloc(N*sizeof(double));
  nv2 = (double*)malloc(N*sizeof(double));
  nt1 = (double*)malloc(N*sizeof(double));
  nt2 = (double*)malloc(N*sizeof(double));
  

  phi = (double*)malloc(N*sizeof(double));
  phiid = (double*)malloc(N*sizeof(double));
  dphi = (double*)malloc(N*sizeof(double));
  dep = (double*)malloc(N*sizeof(double));
  ads = (double*)malloc(N*sizeof(double));
  Vext = (double*)malloc(N*sizeof(double));
  
  param_init();
  
  for (i=0;i<=Rbi;i++) rho[i]=rho0;
  for (i=Rbi+1;i<N;i++) rho[i]=0.0;
 
  for (i=0;i<N;i++) {
    n3[i]=0.0;
    n2[i]=0.0;
    nv2[i]=0.0;
    Vext[i]=0.0;
  }

  for (i=0;i<=2*R;i++) {
    z=(i-R)*dr;
    lw3[i]=dr*pi*(lR*lR-z*z);
    lw2[i]=dr*pi_2*lR;
    lwv2[i]=dr*pi_2*z;
    lwt21[i]=dr*pi_2*0.25*(z*z-lR*lR)*(z*z-lR*lR)/lR;
    lwt22[i]=-dr*pi_2*z*(z*z-lR*lR)/lR;
    lwt23[i]=dr*pi_2*z*z/lR;
  }

  init_matrix(R,w3,lw3);
  init_matrix(R,w2,lw2);
  init_matrix(R,wv2,lwv2);
  init_matrix(R,wt21,lwt21);
  init_matrix(R,wt22,lwt22);
  init_matrix(R,wt23,lwt23);

}

void weighted_densities()
{
  int i,j,jp,jm,k;
  double hr,hrho,n3i,n2i,nv2i,nt2i;

  n3[0]=pi_4*fr2int(rho,dr,R,0.0);
  n2[0]=pi_4*rho[R];
  nv2[0]=0.0;
  nt1[0]=pi_4/3.*rho[R];
  nt2[0]=nt1[0];
  
  for(i=1; i<R; i++) {

    hr = i*dr;
    
    n3i = 0.;
    n2i = 0.;
    nv2i = 0.;
    nt2i = 0.;

    jp = R+i;
    jm = R-i;

    for(j = jm; j <= jp; j++) {
      hrho = rho[j]*j*dr;
      n3i += hrho*w3[jp-jm][R+i-j];
      n2i += hrho*w2[jp-jm][R+i-j];
      nv2i += hrho*wv2[jp-jm][R+i-j];
      nt2i += hrho*(wt21[jp-jm][R+i-j]/(hr*hr*hr)+wt22[jp-jm][R+i-j]/(hr*hr)+wt23[jp-jm][R+i-j]/hr);
    }
   
    n3[i] = pi_4*fr2int(rho,dr,R-i,0.0)+n3i/hr;
    n2[i] = n2i/hr;
    nv2[i] = (n3i/hr+nv2i)/hr;
    nt1[i] = 0.5*(n2i/hr-nt2i);
    nt2[i] = nt2i;
  }

  for(i=R; i<N; i++) {

    hr = i*dr;
    
    n3i = 0.;
    n2i = 0.;
    nv2i = 0.;
    nt2i = 0.;

    jp = i+R;
    jm = i-R;

    if (jp>Rbi) jp=Rbi;

    k=i>=Rbi-R? -1 : 1;

    for(j = jm; j <= jp; j++) {
      hrho = rho[j]*j*dr;
      n3i += hrho*w3[jp-jm][R+k*(i-j)];
      n2i += hrho*w2[jp-jm][R+k*(i-j)];
      nv2i += k*hrho*wv2[jp-jm][R+k*(i-j)];
      nt2i += hrho*(wt21[jp-jm][R+k*(i-j)]/(hr*hr*hr)+k*wt22[jp-jm][R+k*(i-j)]/(hr*hr)+wt23[jp-jm][R+k*(i-j)]/hr);
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
  double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double hR,double a,double b);
  double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double hR,double a,double b);
  double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double hR,double a,double b);
  double mydt1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double hR,double a,double b);
  double mydt2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double hR,double a,double b);
             
  int i,j,jp,jm,k;
  double d2[N],d3[N],dv2[N],dt1[N],dt2[N];
  double lR=R*dr,hr,hrs;
  double c3i,c2i,cv21i,cv22i,ct1,ct2;
  double ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;

  for (i=0;i<N;i++) {
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

  dep[0]=pi_4*fr2int(d3,dr,R,0.0)+pi_4*(d2[R]+dv2[R]+dt2[R])-mu;
  dphi[0]=dep[0]+log(rho[0]/rho0);
  
  for(i=1; i<R; i++) {  

    hr = i*dr;
    c3i = 0.;
    c2i = 0.;
    cv21i = 0.;
    cv22i = 0.;
    ct1 = 0.0;
    ct2 = 0.0;
    
    jp = R+i;
    jm = R-i;

    for(j=jm; j<=jp; j++) {
      
      hrs = j*dr;

      c3i += d3[j]*hrs*w3[jp-jm][R+i-j];
      c2i += d2[j]*hrs*w2[jp-jm][R+i-j];
      cv21i += dv2[j]*w3[jp-jm][R+i-j]; 
      cv22i -= dv2[j]*hrs*wv2[jp-jm][R+i-j];
      ct1 += dt1[j]*hrs*0.5*w2[jp-jm][R+i-j]-dt1[j]*0.5*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
      ct2 += dt2[j]*(wt21[jp-jm][R+i-j]/hrs-wt22[jp-jm][R+i-j]+hrs*wt23[jp-jm][R+i-j]);
    }
    
    dep[i] = pi_4*fr2int(d3,dr,R-i,0.0) + (c3i + c2i + cv22i + cv21i + ct1 + ct2)/hr - mu;
    
    dphi[i] = dep[i]+log(rho[i]/rho0);
    
  } 
 
 for(i=R; i<N; i++) {  

    hr = i*dr;
    c3i = 0.;
    c2i = 0.;
    cv21i = 0.;
    cv22i = 0.;
    ct1 = 0.0;
    ct2 = 0.0;
    
    jp = i+R;
    jm = i-R;

    if (jp>N-1) jp=N-1;

    k = i>=Rbi-1? -1 : 1;

    for(j=jm; j<=jp; j++) {
      
      hrs = j*dr;

      c3i += d3[j]*hrs*w3[jp-jm][R+k*(i-j)];
      c2i += d2[j]*hrs*w2[jp-jm][R+k*(i-j)];

      /*
      w3 part does not have factor of k due to symmetry.
      wv2 however has to be corrected.
      */
      
      cv21i += dv2[j]*w3[jp-jm][R+k*(i-j)]; 
      cv22i -= k*dv2[j]*hrs*wv2[jp-jm][R+k*(i-j)];
    
      if (i==R) ct1 = ct2 = 0.0;

      else {
        ct1 += dt1[j]*hrs*0.5*w2[jp-jm][R+k*(i-j)]-dt1[j]*0.5*(wt21[jp-jm][R+k*(i-j)]/hrs-k*wt22[jp-jm][R+k*(i-j)]+hrs*wt23[jp-jm][R+k*(i-j)]);
        ct2 += dt2[j]*(wt21[jp-jm][R+k*(i-j)]/hrs-k*wt22[jp-jm][R+k*(i-j)]+hrs*wt23[jp-jm][R+k*(i-j)]);
      }
    }
    
    dep[i] = (c3i + c2i + cv22i + cv21i + ct1 + ct2)/hr - mu;
    
    if (i<N-R) dphi[i]=dep[i]+log(rho[i]/rho0);
    else dphi[i]=0.0;
  } 

  //i=R has to be treated carefully since wt21 has otherwise division with zero (hrs=0 due to jm=0 at i=R)
  dep[R] = 0.5*(dep[R-1]+dep[R+1]);
  dphi[R] = 0.5*(dphi[R-1]+dphi[R+1]);

}

double omega(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double hr,ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,sum;
  double lR=R*dr;
  int i;

  weighted_densities();

  for (i=0;i<N;i++) {
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

    if (i<=Rbi) phiid[i]=rho[i]*(log(rho[i]/rho0)-1.0-mu)+mode*p_bulk;
    else phiid[i] = 0.0;
  }

  sum = fr2int(phi,dr,N-1,0.0)+fr2int(phiid,dr,Rbi,0.0);


  if (mode==0) return pi_4*sum;
  else return sum/(Rb*Rb);

}

double surface_tension(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2,double nt1,double nt2,double a,double b);

  double hr,ln3,ln2,ln1,ln0,lnv2,lnv1,lnt1,lnt2;
  double phi1,phi2,phi3,sum;
  double lR=R*dr;
  int i;

  weighted_densities();

  for (i=0;i<N;i++) {
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

    if (rho[i]>0) phiid[i]=rho[i]*(log(rho[i]/rho0)-1.0-mu);
    else phiid[i] = 0.0;
  }

  sum = fr2int(phi,dr,N,0.0)+fr2int(phiid,dr,Rbi,0.0);


  if (mode==0) return pi_4*sum;
  else return sum/(Nr*Nr);

}

double adsorption(int mode) {

   int i;

   for (i=0; i<N; i++) ads[i] = rho[i]-mode*rho0;

   return fr2int(ads,dr,Rbi,0.0)/(Rb*Rb);

}

void minimize(int show)
{
  double hw=0.0,hw_old=1.0,hw_eps,hr;
  int i,j=1;

  hw_eps=1.e-21;

  for (i=0;i<=Rbi;i++) rho[i]=rho0;
  for (i=Rbi+1;i<N;i++) rho[i]=0.0;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<500000)) {

    hw_old=hw;

    weighted_densities();

    derivative();

    for (i=0;i<=Rbi;i++) {
      hr=rho[i];
      rho[i]=(1-a1)*rho[i]+a1*rho0*exp(-dep[i]);
      hw+=(hr-rho[i])*(hr-rho[i]);
    }

    if (show==1) printf("%d %20.6lg\n", j, (hw-hw_old)/hw_eps);

    j++;
  }
}

double adsorption2() {

  double deta=1.e-6,etam=eta,dg=0.0,dmu=0.0;
  int i;


  for (i=0;i<=1;i++) {

    eta = etam-pow(-1,i)*deta;

    param_init();
    minimize(0);

    dmu += pow(-1,i+1)*mu/(2.*deta);
    dg += pow(-1,i+1)*omega(1)/(2.*deta);

  }

  dmu += 1./etam;

  eta = etam;

  return -dg/dmu;

}

double contact() {

  double dgdR=0.0,Rmid=Rb,g=omega(1);
  int i, Rbimid=Rbi;
  
  
  for (i=0;i<=1;i++) {

    Rb = Rmid+pow(-1,i+1)*dr;
    Rbi = Rbimid + pow(-1,i+1);

    minimize(0);
    
    dgdR += pow(-1,i+1)*omega(1)/(2.*dr);

  }

  Rb = Rmid;
  Rbi = Rbimid;

  return p_bulk-2.*g/Rb-dgdR;

}

double vex() {

  for (int i=0; i<N; i++) ads[i] = rho[i]/rho0;

  return Nr/3.-fr2int(ads, dr, Rbi, 0.0)/(Nr*Nr);


}


int main(int argc,char **argv)
{
  void filename(char *fln); 
  int i,j;
  char fln[250],fname[50];
  
   if (argc>=5) {
    Nr=atof(argv[1]);
    eta=atof(argv[2]);
    A=atof(argv[3]);
    B=atof(argv[4]);
    if (argc>=6) a1=atof(argv[5]);
    else a1=0.1;
    
    data_init();
    filename(fname);
   
    minimize(0);
 
    sprintf(fln,"sp_%s_%.2f_%.2f_%.2f_%.2f",fname,1.*N/R,eta,A,B);
    fp=fopen(fln,"w");
    for (i=0;i<N;i++) {
      fprintf(fp,"%15.10lf %20.19e %20.19e\n",i*dr/2.,rho[i],dphi[i]);
    }
    fclose(fp);
    
    /*printf("omega = %.6f\n", omega(0));
    printf("gamma = %.6f\n", 4.*surface_tension(1));
    printf("<N>   = %.6f\n", pi_4*fr2int(rho,dr,Rbi,0.0));
    printf("v_ex  = %.6f\n", 0.5*vex());
    printf("adsorption %.6f, %.6f\n", adsorption(1),adsorption2());
    printf("contact %.6f %.6f\n", rho[Rbi], contact());
    */

    printf("%.4f %lf %lf\n", Nr, 4.*surface_tension(1), 0.5*vex());

    free(rho);
    free(n3);
    free(n2);
    free(nv2);
    free(nt1);
    free(nt2);
 
    free(phi);
    free(phiid);
    free(dphi);
    free(dep);
    free(ads);
    free(Vext);

  } else fprintf(stderr,"syntax: %s Nr eta A B [a]\n",argv[0]);
}


