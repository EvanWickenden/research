
#include <iostream>
#include <assert.h>

#include "lattice2.hpp"

#define SUCCESS()  std::cout << "success!" << std::endl;

#define abs_val(a)  ({ double _a = a; (_a < 0) ? -_a : _a; })


int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    static Lattice2<3,3,3,3> l(gen, 0);

    std::cout << "testing total action" << std::endl;
    assert(2 * 3 * 3 * 3 * 3 * 6 == (int) l.total_action());
    SUCCESS();

    static Lattice2<4,4,4,4> l2(gen, 1);


    int d[] = {1,2,3,4};

    // compilation should fail if below line commented out
//    l2.plaquette<0,0>(d);   

    std::cout << "testing (supposedly) equivalent plaquettes" << std::endl;
    {
        su2 t1 = l2.plaquette<0,1>(d);
        su2 t12 = l2.plaquette<0,1>(d);
        assert(t1 == t12);

        su2 t2 = l2.plaquette<1,0>(d);
        std::cout << t1 << t2 << std::endl;
        assert(t1 == t2);
    }
    SUCCESS();

    std::cout << "testing various actions" << std::endl;

#define test_action(mu, d, precision)                           \
    {                                                           \
        su2 delta = su2::generate2(gen, 1);                     \
        su2& link = l2(d[0], d[1], d[2], d[3])->links[mu];      \
        double local_action = l2.template link_action<mu>(d);            \
        double tot_action = l2.total_action();                  \
        link *= delta;                                          \
        double new_local_action = l2.link_action<mu>(d);        \
        double new_tot_action = l2.total_action();              \
        assert(abs_val(new_local_action - local_action - new_tot_action + tot_action) < precision);\
        SUCCESS();                                              \
    };

    test_action(0, d, 0.0001);
    test_action(1, d, 0.0001);
    test_action(2, d, 0.0001);
    test_action(3, d, 0.0001);

    int d2[] = {0,0,0,0};
    test_action(0, d2, 0.0001);
    test_action(1, d2, 0.0001);
    test_action(2, d2, 0.0001);
    test_action(3, d2, 0.0001);
}
