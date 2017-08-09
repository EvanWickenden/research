
#include <stdexcept>

#include "su2.h"


su2& su2::operator=(const su2& rhs)
{
    a = rhs.a;
    b = rhs.b;
    return *this;
}

su2& su2::operator=(const su2&& rhs)
{
    a = rhs.a;
    b = rhs.b;
    return *this;
}

int su2::operator==(const su2& rhs)
{
    return a == rhs.a && b == rhs.b;
}

su2 su2::inverse() const
{
    return su2(std::conj(a), -b);
}

Complex su2::trace()
{
    return a + std::conj(a);
}

void su2::re_unitarize()
{
    double det = std::norm(a) + std::norm(b);
    a /= det;
    b /= det;
}

su2 su2::identity()
{
    return su2(1,0);   
}

static Complex _i = Complex(0,1);

Matrix<2> su2::u1()
{
    Matrix<2> M =  {{ {0, _i}, {_i, 0} }};
    return M;
}

Matrix<2> su2::u2()
{
    Matrix<2> M = {{ {0, -1}, {1, 0} }};
    return M;
}

Matrix<2> su2::u3()
{
    Matrix<2> M = {{ {_i, 0}, {0, -_i} }};
    return M;
}


#define abs_val(x)  ({ double _x = (x); (_x < 0) ? -_x : _x; })

su2 su2::generate2(std::mt19937& g, double width)
{
    std::normal_distribution<double> _c(0, width);
    Matrix<2> M = exp( _c(g) * su2::u1() + _c(g) * su2::u2() + _c(g) * su2::u3() );
    su2 s(M(0,0), M(1,0));
    s.re_unitarize();
    return s;
}


su2 su2::generate(std::mt19937& g, 
    double av_norm_sq_a, double std_dev, double theta_std_dev)
{
    static std::normal_distribution<double> _theta(0, theta_std_dev);

    std::normal_distribution<double> _nsq_a(av_norm_sq_a, std_dev);

    if (av_norm_sq_a > 1)
    {
        throw std::domain_error("cannot provide average norm squared greater than one");
    }
    if (av_norm_sq_a < 0)
    {
        throw std::domain_error("cannot provide average norm squared less than zero");
    }

    double nsq_a, nsq_b; 

    Complex i_theta_a(0, _theta(g));
    Complex i_theta_b(0, _theta(g));

    while ((nsq_a = _nsq_a(g)) >= 1 || (nsq_a <= 0))
        ;

    nsq_b = 1 - nsq_a;

    double na = sqrt(nsq_a);
    double nb = sqrt(nsq_b);

    return su2(na * std::exp(i_theta_a), nb * std::exp(i_theta_b));
}

std::ostream& operator<<(std::ostream& stream, const su2& matrix)
{
    stream << "/  " << matrix.a << "   " << -std::conj(matrix.b) << "  \\ \n";
    stream << "\\  " << matrix.b << "   " << std::conj(matrix.a)  << "  / \n";
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const su2&& matrix)
{
    stream << "/  " << matrix.a << "   " << -std::conj(matrix.b) << "  \\ \n";
    stream << "\\  " << matrix.b << "   " << std::conj(matrix.a)  << "  / \n";
    return stream;
}


su2& su2::operator *= (const su2& rhs)
{
    Complex _a = a * rhs.a - std::conj(b) * rhs.b;
    Complex _b = b * rhs.a + std::conj(a) * rhs.b;
    a = _a;
    b = _b;
    return *this;
}

su2 operator * (const su2& lhs, const su2& rhs)
{
    return su2(lhs.a, lhs.b) *= rhs;
}

su2& su2::operator *= (const double c)
{
    a *= c;
    b *= c;
    return *this;
}

su2 operator * (const double c, const su2& rhs)
{
    return su2(rhs) *= c;
}



