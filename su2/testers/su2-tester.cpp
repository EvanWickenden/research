
#include "../su2.h"

#include <random>
#include <stdexcept>
#include <iostream>
#include <assert.h>

#define SUCCESS()    std::cout << "success\n" << std::endl

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
    
    su2 matrix(su2::generate(gen));

    double det = std::norm(matrix.a) + std::norm(matrix.b);
    std::cout << "testing determinant; = " << det << std::endl;

    assert(det == 1.0);
    SUCCESS();

    std::cout << "testing multiplication " << std::endl;
    su2 identity(1, 0);
    assert(matrix*identity == matrix);
    SUCCESS();

    std::cout << "testing inverses" << std::endl;
    assert(identity.inverse() == identity);

    su2 a(std::complex<double> (1, 2), std::complex<double> (3, 4));
    su2 a_inverse(std::complex<double> (1, -2), - std::complex<double> (3, 4));

    std::cout << a.inverse() << a_inverse << std::endl;

    assert(a.inverse() == a_inverse);
    SUCCESS();


    std::cout << "expect identity:\n" << su2::generate2(gen, 0.0) << std::endl;
}
