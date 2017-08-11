
#include "../su2.h"

#include <random>
#include <stdexcept>
#include <iostream>
#include <assert.h>

#define SUCCESS()    std::cout << "success\n" << std::endl

#define abs_val(x)  ({ double _x = (x); _x < 0 ? -_x : _x; })

int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    /* invalid range */
    /*
    std::cout << "testing invalid parameters" << std::endl;
    try
    {
        su2::generate(gen, 10, 1);
    }
    catch (std::exception& e)
    {
        std::cout  << e.what() << std::endl;
        goto alive;
    }

    std::cerr << "should have thrown exception" << std::endl;
    exit(1);

alive:
    SUCCESS();
    */
    
    su2 matrix = su2::generate2(gen, 10);
    std::cout << matrix << std::endl;
    double det = std::norm(matrix.a) + std::norm(matrix.b);
    std::cout << "testing determinant; = " << det << std::endl;

    assert( abs_val(det - 1) < 0.0001 );
    SUCCESS();

    std::cout << "testing multiplication " << std::endl;
    su2 identity(1, 0);
    assert(matrix*identity == matrix);
    SUCCESS();

    std::cout << "testing inverses" << std::endl;
    assert(identity.inverse() == identity);
//    assert(matrix * matrix.inverse() == identity);
    std::cout << matrix * matrix.inverse() << std::endl;

    su2 a(std::complex<double> (1, 2), std::complex<double> (3, 4));
    su2 a_inverse(std::complex<double> (1, -2), - std::complex<double> (3, 4));

    std::cout << a.inverse() << a_inverse << std::endl;

    assert(a.inverse() == a_inverse);
    SUCCESS();

    std::cout << "testing trace" << std::endl;
    assert(su2::identity().half_trace() == 1);
    su2 random = su2::generate(gen, 1);
    assert(random.a.real()  == random.half_trace());
    SUCCESS();

    std::cout << "expect identity:\n" << su2::generate2(gen, 0.0) << std::endl;
    

    std::cout << "testing de-unitarization due to rounding errors" << std::endl;
    su2 product = su2::identity();
    for (int i = 0; i < 1000; i++)
    {
        product *= su2::generate2(gen, 1);
    }
    std::cout << product.det() << std::endl;

}

