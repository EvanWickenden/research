
#include <assert.h>

#include "callback/callback.h"
#include "lattice.hpp"
#include "path-integral.hpp"
#include "data-set.h"
#include "observables.hpp"

#define abs_val(a) ({ double _a = a; (_a < 0) ? -_a : _a; })

static const int T = 3, X = 8, Y = 8, Z = 8;

template <int nr_configurations>
struct Snapshot : public Observable<T,X,Y,Z>
{
    char lattices[nr_configurations * sizeof(Lattice<T,X,Y,Z>)]; 
    int index;

#define SIZE  (sizeof (Lattice<T,X,Y,Z>))

    Snapshot() : index(0) {}

    void operator () (const Lattice<T,X,Y,Z>& l)
    {
        memcpy(lattices + (index++) * SIZE, &l, SIZE);
    };

    Lattice<T,X,Y,Z>& get(int index)
    {
        return * (Lattice<T,X,Y,Z> *) (lattices + index * SIZE);
    };

#undef SIZE
};

template <typename _Snapshot> struct WilsonLoop
{
    
};


int main()
{
    double beta = 2;
    double nr_configurations = 100;
    double width = 0.2;

//    Quarks<7> q(nr_configurations);

    std::random_device ran;
    std::mt19937 gen(ran());

    static Lattice<T,X,Y,Z> lattice(gen, 0.1);

    std::cout << sizeof(Lattice<T,X,Y,Z>) << std::endl;

//    path_integral<T,X,Y,Z>(lattice, beta, nr_configurations, width, n, gen);

    /*

    auto pi = [&](double beta, double width)
    {
        AvPlaquette<T,X,Y,Z> a(nr_configurations);
        path_integral<T,X,Y,Z>(lattice, beta, nr_configurations, width, a, gen);
        a.analyze();
        std::cout << a << std::endl;
    };

    pi(1, 0.4);
    pi(2, 0.2);
    pi(3, 0.15);
    pi(4, 0.01);
    pi(5, 0.005);
    pi(6, 0.0005);
    */
}

