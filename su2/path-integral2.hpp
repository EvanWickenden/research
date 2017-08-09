#ifndef __PATH_INTEGRAL_2_HPP__
#define __PATH_INTEGRAL_2_HPP__


#include <random>

#include "lattice/lattice2.hpp"
//#include "path-integral.h"

template <int T, int X, int Y, int Z>
using Observable = Callback<void, const Lattice2<T,X,Y,Z>&>;


template <int size>
struct RandomMatrices
{
    su2 matrices[size];

    RandomMatrices(std::mt19937& generator, double width = 0.1)
    {
        int i;
        for (i = 0; i < size; i += 2)
        {
            matrices[i] = su2::generate2(generator, width);
            matrices[i + 1] = matrices[i].inverse();
        }
    }
    su2& operator()(int i)
    {
        return matrices[i];
    }

    friend std::ostream& operator<<(std::ostream& stream, RandomMatrices& random)
    {
        int i;
        for (i = 0; i < size; i++)
        {
            stream << random.matrices[i] << std::endl;
        }
        return stream;
    }
};


template <int T, int X, int Y, int Z>
void path_integral2(Lattice2<T,X,Y,Z>& lattice, double beta, int nr_configurations, double width, Observable<T,X,Y,Z>& observable, std::mt19937& gen)
{
    const int size = 64;
    RandomMatrices<size> random_matrices(gen, 0.1);
    std::uniform_int_distribution<int> _index(0, size - 1);
    std::uniform_real_distribution<double> _accept(0, 1);

    int k = 0;

    int acceptance_ratio[2] = {0,0};

    for (int n = 0; n < nr_configurations; n++)
    {
        int coords[4];
        su2 copy, delta;
        double old_action, new_action;
        typename Lattice2<T,X,Y,Z>::LatticeSite *site;

        /* old school lambda */
#define propose_change(__n)                                                             \
        {                                                                               \
            su2& link = site->links[__n];                                               \
            delta = random_matrices(_index(gen));                                       \
            copy = link;                                                                \
            old_action = lattice.template link_action<__n>(coords);                     \
            link = delta * link;                                                        \
            new_action = lattice.template link_action<__n>(coords);                     \
            acceptance_ratio[1]++;                                                      \
            if (new_action > old_action                                                 \
                    && _accept(gen) > exp( 0.5 * beta * (old_action - new_action)))     \
                link = copy;                                                            \
            else acceptance_ratio[0]++;                                                 \
        };                  

        int &t = coords[0], &x = coords[1], &y = coords[2], &z = coords[3];

        for (t = 0; t < T; t++) {
            for (x = 0; x < X; x++) {
                for (y = 0; y < Y; y++) {
                    for (z = 0; z < Z; z++) {
                        site = lattice(t,x,y,z);                       
                        propose_change(0);
                        propose_change(1);
                        propose_change(2);
                        propose_change(3);
        }}}}

        observable(lattice);

        if (k++ % 50 == 0) lattice.re_unitarize();
#undef propose_change
    }

    std::cout << "acceptance ratio: " << acceptance_ratio[0] / (double) acceptance_ratio[1] << std::endl;
}


#endif
