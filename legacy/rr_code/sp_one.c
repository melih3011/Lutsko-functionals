#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8192
#define R 250

#define pi 3.14159265358979323846
#define pi_2 6.2831853071795864769   /* 2*pi */
#define pi_4 12.5663706143591729539  /* 4*pi */

double rho[N],n3[N],n2[N],nv2[N];
double phi[N],phiid[N],dphi[N],dep[N],ads[N],Vext[N],hr[N];
double w3[2*R+1][2*R+1],w2[2*R+1][2*R+1],wv2[2*R+1][2*R+1];
double w30[R+1][2*R+1],w20[R+1][2*R+1],wv20[R+1][2*R+1];
double eta,rho0,mu,dr=1.0/(1.0*R),a1,p_bulk,Rb;
int Rbi;

FILE *fp;

double integrate_array_section(double f[], double dx, int i0, int i1)
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

double fr2int2(double f[],int n,double ldr, double a)
{
  double val=0.0,hr1=6*a*a,hr2=4*a*ldr,hr3=ldr*ldr;
  double hf[N];
  int j;

  for (j=0;j<=n;j++) {
    hr1=a+j*ldr;
    hf[j]=hr1*hr1*f[j];
  }
    
  return integrate_array_section(hf,ldr,0,n);
}

double fr2int(double f[],int n,double ldr, double a)
{
  double val=0.0,hr1=6*a*a,hr2=4*a*ldr,hr3=ldr*ldr;
  register int j;

  for (j=0;j<n;j++)
    val+=f[j]*(hr1+hr2*(1+3*j)+hr3*(1+4*j+6*j*j))+
      f[j+1]*(hr1+hr2*(2+3*j)+hr3*(3+8*j+6*j*j));
  return(ldr*val/12.0);
}

void data_init()
{
  double pressure(double eta,double rho0);
  double chemical_potential(double eta);
  double lw3[2*R+1],lw2[2*R+1],lwv2[2*R+1];
  int i,j;
  double z,lR=R*dr;

  rho0=3.0*eta/(pi_4);

  mu=chemical_potential(eta);

  p_bulk=pressure(eta,rho0);

  for (i=0;i<=Rbi-R;i++) rho[i]=rho0;
  for (i=Rbi-R+1;i<N;i++) rho[i]=0.0;

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
  }

  for (i=0;i<5;i++)
    for (j=0;j<2*R;j++) {
      w3[i][j]=0.0;
      w2[i][j]=0.0;
      wv2[i][j]=0.0;
    }
  
  w3[1][0]=0.5*lw3[0]; w3[1][1]=0.5*lw3[1];
  w3[2][0]=1./3.*lw3[0]; w3[2][1]=4./3.*lw3[1]; w3[2][2]=1./3.*lw3[2];
  w3[3][0]=3./8.*lw3[0]; w3[3][1]=9./8.*lw3[1]; w3[3][2]=9./8.*lw3[2];
  w3[3][3]=3./8.*lw3[3];
  w3[4][0]=14./45.*lw3[0]; w3[4][1]=64./45.*lw3[1]; w3[4][2]=24./45.*lw3[2];
  w3[4][3]=64./45.*lw3[3]; w3[4][4]=14./45.*lw3[4];

  w2[1][0]=0.5*lw2[0]; w2[1][1]=0.5*lw2[1];
  w2[2][0]=1./3.*lw2[0]; w2[2][1]=4./3.*lw2[1]; w2[2][2]=1./3.*lw2[2];
  w2[3][0]=3./8.*lw2[0]; w2[3][1]=9./8.*lw2[1]; w2[3][2]=9./8.*lw2[2];
  w2[3][3]=3./8.*lw2[3];
  w2[4][0]=14./45.*lw2[0]; w2[4][1]=64./45.*lw2[1]; w2[4][2]=24./45.*lw2[2];
  w2[4][3]=64./45.*lw2[3]; w2[4][4]=14./45.*lw2[4];

  wv2[1][0]=0.5*lwv2[0]; wv2[1][1]=0.5*lwv2[1];
  wv2[2][0]=1./3.*lwv2[0]; wv2[2][1]=4./3.*lwv2[1]; wv2[2][2]=1./3.*lwv2[2];
  wv2[3][0]=3./8.*lwv2[0]; wv2[3][1]=9./8.*lwv2[1]; wv2[3][2]=9./8.*lwv2[2];
  wv2[3][3]=3./8.*lwv2[3];
  wv2[4][0]=14./45.*lwv2[0];wv2[4][1]=64./45.*lwv2[1];wv2[4][2]=24./45.*lwv2[2];
  wv2[4][3]=64./45.*lwv2[3]; wv2[4][4]=14./45.*lwv2[4];
  
  for (i=5;i<=2*R;i++) {
    for (j=0;j<2*R-i;j++) {
      w3[i][2*R-j]=0.0;
      w2[i][2*R-j]=0.0;
      wv2[i][2*R-j]=0.0;
    }
      w3[i][i]=3./8.*lw3[i];
      w3[i][i-1]=7./6.*lw3[i-1];
      w3[i][i-2]=23./24.*lw3[i-2];
      w2[i][i]=3./8.*lw2[i];
      w2[i][i-1]=7./6.*lw2[i-1];
      w2[i][i-2]=23./24.*lw2[i-2];
      wv2[i][i]=3./8.*lwv2[i];
      wv2[i][i-1]=7./6.*lwv2[i-1];
      wv2[i][i-2]=23./24.*lwv2[i-2];
      for (j=3;j<=i-3;j++) {
	w3[i][j]=lw3[j];
	w2[i][j]=lw2[j];
	wv2[i][j]=lwv2[j];
      }
      w3[i][2]=23./24.*lw3[2];
      w3[i][1]=7./6.*lw3[1];
      w3[i][0]=3./8.*lw3[0];
      w2[i][2]=23./24.*lw2[2];
      w2[i][1]=7./6.*lw2[1];
      w2[i][0]=3./8.*lw2[0];
      wv2[i][2]=23./24.*lwv2[2];
      wv2[i][1]=7./6.*lwv2[1];
      wv2[i][0]=3./8.*lwv2[0];
  }

  for (i=0;i<=R;i++)
    for (j=0;j<=2*R+1;j++) {
      w20[i][j]=0.;
      w20[i][j]=0.;
      wv20[i][j]=0.;
    }
  
  for (i=0;i<=R;i++) {
    for (j=R-i;j<=R+i;j++) {
	w30[i][R-i+j]=lw3[R-i+j];
	w20[i][R-i+j]=lw2[R-i+j];
	wv20[i][R-i+j]=lwv2[R-i+j];
    }
  }
  w20[1][2*R]*=1./3.; w20[1][2*R-1]*=4./3.; w20[1][2*R-2]*=1./3.;
  w20[2][2*R]*=14./45.;w20[2][2*R-1]*=64./45.;w20[2][2*R-2]*=24./45.;
  w20[2][2*R-3]*=64./45.; w20[2][2*R-4]*=14./45.;

  w30[1][2*R]*=1./3.; w30[1][2*R-1]*=4./3.; w30[1][2*R-2]*=1./3.;
  w30[2][2*R]*=14./45.;w30[2][2*R-1]*=64./45.;w30[2][2*R-2]*=24./45.;
  w30[2][2*R-3]*=64./45.; w30[2][2*R-4]*=14./45.;

  wv20[1][2*R]*=1./3.; wv20[1][2*R-1]*=4./3.; wv20[1][2*R-2]*=1./3.;
  wv20[2][2*R]*=14./45.;wv20[2][2*R-1]*=64./45.;wv20[2][2*R-2]*=24./45.;
  wv20[2][2*R-3]*=64./45.; wv20[2][2*R-4]*=14./45.;

  
  for (i=3;i<=R;i++) {
    w20[i][2*R]*=3./8.;w20[i][2*R-1]*=7./6.;w20[i][2*R-2]*=23./24.;
    w20[i][2*(R-i)]*=3./8.;w20[i][2*(R-i)+1]*=7./6.;w20[i][2*(R-i)+2]*=23./24.;
    w30[i][2*R]*=3./8.;w30[i][2*R-1]*=7./6.;w30[i][2*R-2]*=23./24.;
    w30[i][2*(R-i)]*=3./8.;w30[i][2*(R-i)+1]*=7./6.;w30[i][2*(R-i)+2]*=23./24.;
    wv20[i][2*R]*=3./8.;wv20[i][2*R-1]*=7./6.;wv20[i][2*R-2]*=23./24.;
    wv20[i][2*(R-i)]*=3./8.;wv20[i][2*(R-i)+1]*=7./6.;wv20[i][2*(R-i)+2]*=23./24.;
  }
  
}

void weighted_densities()
{
  int i, j, jp, jm;
  double hr,hrho,n3i,n3i0,n2i,nv2i,lR=R*dr;
  double eps=1e-12;

  n2[0] = pi_4*lR*lR*rho[R];
  n3[0] = pi_4*fr2int(rho,R,dr,0);
  nv2[0] = 0.0;

  for(i = 1; i <=R ; i++) {
    hr = i*dr;
    
    n3i0 = pi_4*fr2int(rho,R-i,dr,0);
    n3i = 0.;
    n2i = 0.;
    nv2i = 0.;

    for(j = R-i; j <= i+R; j++) {
      hrho = rho[j]*j*dr;
      n3i += hrho*w30[i][R-i+j];
      n2i += hrho*w20[i][R-i+j];
      nv2i += hrho*wv20[i][R-i+j];
    }
    
    n3[i] = n3i0+n3i/hr;
    n2[i] = n2i/hr;
    nv2[i] = (n3i/hr-nv2i)/hr;
  }
  
  for(i = R+1; i < Rbi ; i++) {

    hr = i*dr;
    
    n3i = 0.;
    n2i = 0.;
    nv2i = 0.;

    jp=i+R;
    jm=i-R;
    if (jp>=Rbi-R) jp=Rbi-R;
    
    for(j = jm; j <= jp; j++) {
      hrho = rho[j]*j*dr;
      n3i += hrho*w3[jp-jm][R-i+j];
      n2i += hrho*w2[jp-jm][R-i+j];
      nv2i += hrho*wv2[jp-jm][R-i+j];
    }
    
    n3[i] = n3i/hr;
    n2[i] = n2i/hr;
    nv2[i] = (n3[i]-nv2i)/hr;
    if (fabs(nv2[i])<eps) nv2[i]=0.;
  }
  
  n3[Rbi]= 0.0;
  n2[Rbi]= 0.0;
  nv2[Rbi]= 0.0;
}

void derivative()
{
  double myd2(double n0,double n1,double n2,double n3,double nv1,double nv2,
	      double hR);
  double myd3(double n0,double n1,double n2,double n3,double nv1,double nv2,
	      double hR);
  double mydv2(double n0,double n1,double n2,double n3,double nv1,double nv2,
	       double hR);
  int i,j,ja,jp,jm;
  double d2[N],d3[N],dv2[N];
  double lR=R*dr,hr,hrs;
  double c3i,c3i0,c2i,cv2i,cv2i0;
  double ln3,ln2,ln1,ln0,lnv2,lnv1;
  double eps=1e-12;

  for (i=0;i<=Rbi;i++) {
    ln3=n3[i];
    ln2=n2[i];      
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);  

    d2[i]=myd2(ln0,ln1,ln2,ln3,lnv1,lnv2,lR);
    d3[i]=myd3(ln0,ln1,ln2,ln3,lnv1,lnv2,lR);
    dv2[i]=mydv2(ln0,ln1,ln2,ln3,lnv1,lnv2,lR);
  }
  for (i=Rbi+1;i<N;i++) {
    d3[i]=0.0;
    d2[i]=0.0;
    dv2[i]=0.0;
  }


  c2i = pi_4*lR*lR*d2[R];
  c3i = pi_4*fr2int(d3,R,dr,0);
  cv2i = -pi_4*lR*lR*dv2[R];
  
  dep[0] = (c2i + c3i + cv2i) - mu;

  for(i = 1; i <=R ; i++) {
    hr = i*dr;
    
    c3i0 = pi_4*fr2int(d3,R-i,dr,0);
    cv2i0 = -pi_4*(lR-hr)*(lR-hr)*dv2[R-i];
    c3i = 0.;
    c2i = 0.;
    cv2i = 0.;

    for(j = R-i; j <= i+R; j++) {
      hrs = j*dr;
      c3i += d3[j]*hrs*w30[i][R-i+j];
      c2i += d2[j]*hrs*w20[i][R-i+j];
      cv2i += dv2[j] * (hrs* wv2[2*i][R-j+i] - w3[2*i][R-j+i]);
      //cv2i += dv2[j]*(hrs*wv20[i][R+i-j]-w30[i][R+i-j]);
      //      cv2i += dv2[j]*(-w30[i][R-i+j]+hrs*wv20[i][3*R-i-j]);
    }

    dep[i] = c3i0 + cv2i0*0 + (c3i + c2i + cv2i)/hr - mu;
  }

  for(i = R+1; i < Rbi; i++) {  

    hr = i*dr;
    c3i = 0.;
    c2i = 0.;
    cv2i = 0.;
    
    jp = i + R;
    jm = i - R;
    if (jp>Rbi) jp=Rbi;
    
    for(j = jm; j <= jp; j++) {
      hrs = j*dr;

      c3i += d3[j] * hrs * w3[jp-jm][R-j+i];
      c2i += d2[j] * hrs * w2[jp-jm][R-j+i];
      cv2i += dv2[j] * (hrs* wv2[jp-jm][R-j+i] - w3[jp-jm][R-j+i]);
    }

    dep[i] = (c3i + c2i + cv2i)/hr - mu;
    
    if (rho[i]>0) dphi[i]=dep[i]+log(rho[i]/rho0);
    else dphi[i]=0.0;
  } 
}

double omega(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2);

  double in3,in2,inv2,hr,ln3,ln2,ln1,ln0,lnv2,lnv1;
  double phi1,phi2,phi3,sum;
  double lR=R*dr;
  int i;
  
  weighted_densities();

  for (i=0;i<=Rbi;i++) {
    ln3=n3[i];
    ln2=n2[i];
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);

    phi1=myphi1(ln0,ln1,ln2,ln3,lnv1,lnv2);
    phi2=myphi2(ln0,ln1,ln2,ln3,lnv1,lnv2);
    phi3=myphi3(ln0,ln1,ln2,ln3,lnv1,lnv2);

    phi[i]=(phi1+phi2+phi3)+mode*p_bulk;

    if (rho[i]>0) phiid[i]=rho[i]*(log(rho[i]/rho0)-1.0-mu);
    else phiid[i]=0.0;
  }

  sum = fr2int(phi,Rbi,dr,0.)+fr2int(phiid,Rbi-R,dr,0.);

  if (mode==0) return pi_4*sum;
  else return sum/(Rb*Rb);
}

double omega2(int mode)
{
  double myphi1(double n0,double n1,double n2,double n3,double nv1,double nv2);
  double myphi2(double n0,double n1,double n2,double n3,double nv1,double nv2);
  double myphi3(double n0,double n1,double n2,double n3,double nv1,double nv2);

  double in3,in2,inv2,hr,ln3,ln2,ln1,ln0,lnv2,lnv1;
  double phi1,phi2,phi3,sum;
  double lR=R*dr;
  int i;
  
  weighted_densities();

  for (i=0;i<=Rbi;i++) {
    ln3=n3[i];
    ln2=n2[i];
    ln1=ln2/(pi_4*lR);
    ln0=ln2/(pi_4*lR*lR);
    lnv2=nv2[i];
    lnv1=lnv2/(pi_4*lR);

    phi1=myphi1(ln0,ln1,ln2,ln3,lnv1,lnv2);
    phi2=myphi2(ln0,ln1,ln2,ln3,lnv1,lnv2);
    phi3=myphi3(ln0,ln1,ln2,ln3,lnv1,lnv2);

    phi[i]=(phi1+phi2+phi3)+mode*p_bulk;

    if (rho[i]>0) phiid[i]=rho[i]*(log(rho[i]/rho0)-1.0-mu);
    else phiid[i]=0.0;
  }

  sum = fr2int2(phi,Rbi,dr,0.)+fr2int2(phiid,Rbi-R,dr,0.);

  if (mode==0) return pi_4*sum;
  else return sum/(Rb*Rb);
}

double adsorption(int mode)
{
  int i;
  double sum,lR=R*dr;

  for (i=0;i<=Rbi-R;i++) ads[i]=rho[i]-mode*rho0;

  sum = fr2int(ads,Rbi-R,dr,0.0) - mode*rho0*lR*(lR*lR-3.*lR*Rb+3*Rb*Rb)/3.;

  if (mode==0) return pi_4*sum;
  else return sum/(Rb*Rb);
}

double adsorption2(int mode)
{
  int i;
  double sum,lR=R*dr;

  for (i=0;i<=Rbi-R;i++) ads[i]=rho[i]-mode*rho0;

  sum = fr2int2(ads,Rbi-R,dr,0.0) - mode*rho0*lR*(lR*lR-3.*lR*Rb+3*Rb*Rb)/3.;

  if (mode==0) return pi_4*sum;
  else return sum/(Rb*Rb);
}

double vex()
{
  int i;

  for (i=0;i<=Rbi-R;i++) ads[i]=rho[i]/rho0;

  return Rb/3.-fr2int(ads,Rbi-R,dr,0.0)/(Rb*Rb);
}

void minimize()
{
  double da1,hw=0.0,hw_old=1.0,hw_eps,hr;
  int i,j=1;

  hw_eps=0.000000000000000000001;

  while ((fabs(hw-hw_old)>hw_eps)&&(j<10000)) {
  //  for (j=0;j<500;j++) {
    hw_old=hw;

    // hw=omega(0);

    weighted_densities();

    derivative();
    hw = 0.0;
    for (i=0;i<=Rbi-R;i++) {
      hr=rho[i];
      rho[i]=(1-a1)*rho[i]+a1*rho0*exp(-dep[i]);
      hw+=(hr-rho[i])*(hr-rho[i]);
    }
    //    printf("%i  %lg  %lg\n",j,hw,(hw-hw_old)/hw_eps);
    j++;
  }
}

int main(int argc,char *argv[])
{
  void filename(char *fln);
  
  int i,j;
  char fln[250],fname[50];
  double r,hw;

  if (argc>=3) {
    Rb=atof(argv[1]);
    eta=atof(argv[2]);
    if (argc>=4) a1=atof(argv[3]);
    else a1=0.1;

    Rbi=(int)(R*Rb);
    data_init();
    //weighted_densities();
    //derivative();
    minimize();

    filename(fname);
    
    sprintf(fln,"sp_%s_%10.8lf_%6.4lf",fname,eta,Rb);
    fp=fopen(fln,"w");
    fprintf(fp,"# Rb      = %10.8lf\n",Rb);
    fprintf(fp,"# p_bulk  = %20.19e\n",p_bulk);
    fprintf(fp,"# mu      = %20.19e\n",mu + log(eta));
    fprintf(fp,"# sigma   = %20.19e\n",omega(1));
    fprintf(fp,"# Omega   = %20.19e\n",omega(0));
    fprintf(fp,"# <N>     = %20.19e\n",adsorption(0));
    fprintf(fp,"# Gamma   = %20.19e\n",adsorption(1));
    fprintf(fp,"# vex     = %20.19e\n",vex());
    for (i=0;i<=Rbi;i++) {
      fprintf(fp,"%15.10lf   %20.19e  %20.19e %20.19e %20.19e %20.19lf %20.19lf\n",i*dr,rho[i]/rho0,n3[i],n2[i],nv2[i],dep[i],dphi[i]);
    }
    fclose(fp);

    printf("%6.4lf  %20.19lf  %20.19lf\n",Rb/2.0,4.*omega(1),vex()/2.0);
    
    /*
    fp=fopen("w0.fil","w");
    for (i=0;i<=R;i++)
      for (j=0;j<=2*R;j++)
	fprintf(fp,"%i %i %e %e %e\n",i,j,w20[i][j],w30[i][j],wv20[i][j]);
    fclose(fp);
    */
  } else fprintf(stderr,"syntax: %s Rb eta [a]\n",argv[0]);
}


