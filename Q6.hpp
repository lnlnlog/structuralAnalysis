/*
 * Q6.hpp
 * delcare local bond order parameter calculation function, w()
 */

#ifndef Q6_HPP
#define Q6_HPP

#include <complex>

using std::complex;
using std::norm;

const int MAXBOND = 100; // maximum number of bonds to be processed

class Q6 {
private:
    static bool fzero(const double); // test if a double number is close to zero
    static complex<double> sphHarmonic(int, int, double, double); // spherical harmonic functions
    static long factorial(int); // factorial function
    static double legendre(int, int, double); // legendre function

public:
    static int w(const int, const double*, const double*, const double*, const double*, const int, complex<double> [][13], const int); // w function only do q6 not wj
};

#endif
