
#include <iostream>
#include <assert.h>

#include "../lattice.hpp"
#include "../su2.h"

#define SUCCESS()  std::cout << "success!" << std::endl;

#define abs_val(a)  ({ double _a = a; (_a < 0) ? -_a : _a; })


int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    static Lattice<3,3,3,3> l(gen, 0);

    std::cout << "testing total action" << std::endl;
    assert(- 3 * 3 * 3 * 3 * 6 == (int) l.total_action());
    SUCCESS();

    static Lattice<4,4,4,4> l2(gen, 1);


    int d[] = {1,2,3,4};

    // compilation should fail if below line commented out
//    l2.plaquette<0,0>(d);   

    std::cout << "testing plane-inverted plaquettes" << std::endl;
    {
        su2 t1 = l2.plaquette<0,1>(d);
        su2 t12 = l2.plaquette<1,0>(d);
        assert(approx_equal(t1, t12.inverse()));
        std::cout << t1 << t12 << std::endl;
    }
    SUCCESS();

    std::cout << "testing local actions about various links" << std::endl;
    std::uniform_int_distribution<int> _index(0, 3);

#define test_action(mu, d, precision)                           \
    {                                                           \
        d[0] = _index(gen); d[1] = _index(gen); d[2] = _index(gen); d[3]  = _index(gen); \
        su2 delta = su2::generate2(gen, 1);                     \
        su2& link = l2(d[0], d[1], d[2], d[3])->links[mu];      \
        double local_action = l2.template link_action<mu>(d);   \
        double tot_action = l2.total_action();                  \
        link *= delta;                                          \
        double new_local_action = l2.link_action<mu>(d);        \
        double new_tot_action = l2.total_action();              \
        assert(abs_val(new_local_action - local_action - new_tot_action + tot_action) < precision);\
    };

    test_action(0, d, 0.0001);
    test_action(1, d, 0.0001);
    test_action(2, d, 0.0001);
    test_action(3, d, 0.0001);
    SUCCESS();                                              

    std::cout << "testing willson loops" << std::endl;

#define test_wilson(mu, nu, d)                                                              \
    {                                                                                       \
        d[0] = _index(gen); d[1] = _index(gen); d[2] = _index(gen); d[3]  = _index(gen);    \
        assert( approx_equal(l2.wilson_loop<mu, nu>(1,1,d), l2.plaquette<mu, nu>(d)) );     \
    }

    test_wilson(0,1, d);
    test_wilson(0,2, d);
    test_wilson(0,3, d);
    test_wilson(1,2, d);
    test_wilson(1,3, d);
    test_wilson(2,3, d);
    SUCCESS();

    {
#define factor 7.3
        struct FalseLattice : Callback<su2>
        {
            su2 operator () () 
            {
                return factor * su2::identity();
            }
        };
        Lattice<4,4,4,4> l((FalseLattice()));

        int length = 4, width = 3;

        assert(approx_equal(pow(factor, 2*(length + width)), l.wilson_loop<0,1>(length, width, d).half_trace()));
#undef factor
    }
    SUCCESS();


    st::cout << "testing piecewise Wilson Loop vs pre-contracted Wilson Loop" << std::endl;


}



