/* 
 * Q6.cpp
 * implementation of w function
 */

#include <cmath>
#include <iostream>
#include "const.hpp"
#include "Q6.hpp"
#include <cstdlib>

using std::cout;
using std::endl;

int Q6::w(const int l, const double* px, const double* py, const double* pz, const double* area, const int nnp, complex<double> qcmplx[][13], const int partTag) { // w function only do q6 not wj
    
    if(nnp <= 0) return 0;
    
    double cosTheta[MAXBOND], phi[MAXBOND];
    
    if(nnp > MAXBOND) {
        cout << "Error in w(): number of bonds exceeds MAXBOND" << endl;
        exit(1);
    }
    
    for(int i=0; i<=nnp-1; i++) {
        if(fzero(px[i])) {
            if(fzero(py[i])) phi[i] = 0.0;
            else if(py[i] > 0) phi[i] = PI / 2.0;
            else phi[i] = 3.0 / 2.0 * PI;
        }
        else {
            phi[i] = atan(py[i] / px[i]);
            if(px[i] < 0.0) phi[i] += PI;
            else if(py[i] < 0.0) phi[i] += 2.0 * PI;
        }
    
        cosTheta[i] = pz[i] / sqrt( px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i] );
    }
    
    complex<double> Q[2 * l + 1];
    complex<double> c;
    double totalArea = 0.0;
    
    for(int j=0; j<=nnp-1; j++) {
        Q[l] += sphHarmonic(l, 0, cosTheta[j], phi[j]) * area[j]; // m = 0
        for(int m=1; m<=l; m++) {
            c = sphHarmonic(l, m, cosTheta[j], phi[j]) * area[j];
            Q[m+l] += c;
            double sgn;
            m % 2 == 0 ? sgn = 1.0 : sgn = -1.0;
            Q[-m + l] += complex<double> ( sgn * c.real(), -sgn * c.imag());
        }
        totalArea += area[j];
    }
    
    for(int m=-l; m<=l; m++) Q[m + l] /= totalArea; // Q[m + l] /= static_cast<double> (nnp);
    for(int m=-l; m<=l; m++) qcmplx[partTag][m + l] = Q[m + l]; // store all qlm into 2d complex array

    return 0;
}

bool Q6::fzero(const double a) { // test if a double number is close to zero
    return fabs(a) < 1.0e-30;
}

complex<double> Q6::sphHarmonic(int l, int m, double cosTheta, double phi) { // spherical harmonic functions
    int m1 = abs(m);
    
    double c;

    c = sqrt( ( 2.0 * l + 1.0 ) * factorial( l-m1 ) / ( 4.0 * PI * factorial( l+m1 ) ) );
    
    c *= legendre(l, m1, cosTheta);
    
    complex<double> y;
    
    if(m<0) {
        y = complex<double> (cos(m1 * phi), -sin(m1 * phi)) * ((m&1) ? -1.0 : 1.0);
    }
    else {
        y = complex<double> (cos(m1 * phi), sin(m1 * phi));
    }
    
    y *= c;
    
    return y;
}

long Q6::factorial(int n) { // factorial function
    if( n==0 ) return 1;
    else if( n==1 ) return 1;
    else if( n==2 ) return 2;
    else if( n==3 ) return 6;
    else if( n==4 ) return 24;
    else if( n==5 ) return 120;
    else if( n==6 ) return 720;
    else if( n<0 )
    {
        long f=1;
        for( long i=-1; i>=n; i-- ) f *= i;
        return f;
    }
    else
    {
        long f = 720;
        for( long i=7; i<=n; i++ ) f *= i;
        return f;
    }
}

double Q6::legendre(int l, int m, double x) { // legendre function
    if(m < 0 || m > l || fabs(x) > 1.0) {
        cout << "Error in legendre: bad argument" << endl;
        exit(1);
    }
    
    double pmm,somx2,fact;
    
    pmm = 1.0;
    
    if( m > 0 ) {
        somx2 = sqrt( ( 1.0 - x ) * ( 1.0 + x ) );
        fact = 1.0;
        for(int i=1; i<=m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    
    double pmmp1,pll;
    
    if(l == m) return pmm;
    else {
        pmmp1 = x * ( 2*m + 1 ) * pmm;
        if( l == ( m + 1 ) ) return pmmp1;
        else {
            for(int ll=m+2; ll<=l; ll++) {
                pll = ( x * ( 2*ll -1 ) * pmmp1 - ( ll + m - 1 ) * pmm ) / ( ll - m );
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}
