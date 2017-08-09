
#include <iostream>
#include <assert.h>

#include "../path-integral.h"

#define SUCCESS()   std::cout << "success!" << std::endl;

int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    std::normal_distribution<double> c(0, 0.1);

    Matrix<2> M = exp( c(gen)*su2::u1() + c(gen)*su2::u2() + c(gen)*su2::u3() ); 

    std::cout << "testing exponentially generated su2 matrix" << std::endl;
    assert(M(0,0) == std::conj(M(1,1)));
    assert(M(1,0) == -std::conj(M(0,1)));
    SUCCESS();
    
    RandomMatrices random(gen);

    std::cout << "randomly generated matrices" << std::endl;
    std::cout << random << std::endl;
}
