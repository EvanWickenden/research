
#include <iostream>
#include <assert.h>

#include "matrix.h"
#include "tuple.h"

#define SUCCESS()    std::cout << "success!" << std::endl;

int main()
{
    Matrix<2> M = {{ {0, 1}, {1, 0} }};

    Matrix<2> I = Matrix<2>::identity();

    std::cout << "testing equality" << std::endl;
    assert( M == M );
    SUCCESS();

    std::cout << "testing identity element" << std::endl;
    assert( M == M * I );
    assert( I * M == M );
    SUCCESS();

    std::cout << "testing multiplication" << std::endl;
    assert( I == M * M );
    assert( M == M * M * M );
    SUCCESS();

    Matrix<2> A = {{ {1, 2}, {2, 1} }};
    Matrix<2> B = {{ {0, -2}, {-2, 0} }};

    std::cout << "testing addition" << std::endl;
    assert( I == A + B );
    assert( B + A == A + B );
    SUCCESS();

    std::cout << "testing exp(I)" << std::endl;
    std::cout << exp(I) << std::endl;
    std::cout << exp(I, 20) << std::endl;

    Matrix<2> zero;
    assert(I == exp(zero));
    SUCCESS();
}
