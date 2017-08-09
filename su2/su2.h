
#ifndef __SU2_H__
#define __SU2_H__

#include <complex>
#include <iostream>
#include <random>
#include <math.h>

#include "matrix/matrix.h"

typedef std::complex<double> Complex;


/**
 *   / a  -b* \
 *   \ b   a* /
 */

struct su2
{
    Complex a, b;

    static su2 identity();
    static Matrix<2> u1();
    static Matrix<2> u2();
    static Matrix<2> u3();

    static su2 generate(std::mt19937& g, double av_norm_sq_a = 0.5, double std_dev = 0.5, double theta_std_dev = 2 * M_PI);
    static su2 generate2(std::mt19937& g, double width = 0.01); 

    su2() {}
    su2(Complex a, Complex b) : a(a), b(b) {}
    su2(const su2& rhs) : a(rhs.a), b(rhs.b) {}
    su2(const su2&& rhs) : a(rhs.a), b(rhs.b) {}
    su2& operator=(const su2& rhs);
    su2& operator=(const su2&& rhs);

    su2 inverse() const;
    void re_unitarize();


    su2& operator *= (const su2& rhs);
    su2& operator *= (const double c);
    int operator == (const su2& rhs);
    
    Complex trace();

    friend std::ostream& operator<<(std::ostream& stream, const su2& matrix);
    friend std::ostream& operator<<(std::ostream& stream, const su2&& matrix);
};

su2 operator * (const su2& lhs, const su2& rhs);
su2 operator * (const double c, const su2& rhs);

#endif 
