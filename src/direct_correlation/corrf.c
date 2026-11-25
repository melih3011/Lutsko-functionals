#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dx 0.01

double eta, A, B;


double pyCorr(double e, int i) {
    
    double x = i*dx;

    return -pow(1.+2.*e,2)/pow(1.-e,4) + 3.*e*pow(2.+e,2)/(2.*pow(1.-e,4))*x - e*pow(1.+2.*e,2)/(2.*pow(1.-e,4))*pow(x,3);

}


double lkCorr(double a, double b, double e, int i) {
    
    double first, second;

    double x = i*dx;

    first = e/(2.*pow(1.-e,2))*x*(2.*(a+b)*x*x+3.-2.*b-4.*a);
    second = -(8.*a+2.*b-9.)*e*e/pow(1.-e,3)*(1.-x)*(1.+0.5*e/(1.-e)*(1.-x)*(2.+x));

    return pyCorr(e, i) + first + second;

}

int main(int argc, char** argv) {
    
    eta = atof(argv[1]);
    
    for (int j=0; j<=(int)(1./dx); j++) {
        
        printf("%lf %lf %lf %lf %lf\n", j*dx, pyCorr(eta, j), lkCorr(1.3, -1.0, eta, j), lkCorr(1.0, 0.0, eta, j), lkCorr(0.0, 4.0, eta, j));

    }

    return 0;
}
