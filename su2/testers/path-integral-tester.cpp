
#include <iostream>
#include <assert.h>

#include "../path-integral.hpp"
#include "../lattice.hpp"
#include "../su2.h"


int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    
    const int size = 100;
    static RandomMatrices<size> random_matrices(gen, 0.01);

    std::cout << random_matrices << std::endl;

    /*
    std::cout << "testing inverses" << std::endl;
    for (int i = 0; i < size; i+=2)
    {
        su2 product = random_matrices.matrices[i] * random_matrices.matrices[i + 1];
        std::cout << product;
    }
    */
}

