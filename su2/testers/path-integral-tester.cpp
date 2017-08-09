
#include <math.h>
#include <assert.h>

#include "../path-integral.h"


#define SUCCESS()   std::cout << "success!" << std::endl;



int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    int i;
    int N = 10000;

    std::uniform_real_distribution<double> _accept(0, 1);

    Lattice lattice((int[]) {3,3,3,3}, gen);
    LinkAction action = &Lattice::t_link_action;

    su2& link = lattice(0,0,0,0)->forward_t;
    su2 delta = su2::generate(gen);

    std::cout << "testing change in action evalutation" << std::endl;

    double old_tot_action, new_tot_action;
    double old_action, new_action;

    old_tot_action = lattice.total_action();
    old_action = (lattice.*action)(0,0,0,0);
    link *= delta;
    new_tot_action = lattice.total_action();
    new_action = (lattice.*action)(0,0,0,0);

    std::cout << "delta total lattice action " << new_tot_action - old_tot_action << std::endl;
    std::cout << "delta local action " << new_action - old_action << std::endl;

    double beta = 2;

    std::cout << "testing update rule" << std::endl;

    su2 copy(link);

    delta = su2::generate(gen);
    old_action = (lattice.*action)(0,0,0,0);
    link *= delta;
    new_action = (lattice.*action)(0,0,0,0);
    link = copy;

    double probability = exp( 0.5 * beta * (old_action - new_action));

    std::cout << "old action: " << old_action << "; new action: " << new_action << std::endl;
    std::cout << "probability: " << probability << std::endl;

    ProposeChange change(lattice, 2, 0,0,0,0);

    for (i = 0; i < N; i++)
    {
        double accept = _accept(gen);
        change(link, action, delta, accept);
        link = copy;
    }

    std::cout << "acceptance_ratio: " << lattice.acceptance_ratio[0] / (double) lattice.acceptance_ratio[1] << "\n" << std::endl;


    std::cout << "testing change in local vs total action for each link direction" << std::endl;

    std::uniform_int_distribution<int> _index(0,2);

    std::cout << "t" << std::endl;
    {
        int i = _index(gen);
        int j = _index(gen);
        int k = _index(gen);
        int l = _index(gen);

        su2 delta = su2::generate(gen);
        su2& link = lattice(i,j,k,l)->forward_t;

        double old_action, new_action;
        double old_total, new_total;

        old_action = lattice.t_link_action(i,j,k,l);
        old_total = lattice.total_action();

        link *= delta;

        new_action = lattice.t_link_action(i,j,k,l);
        new_total = lattice.total_action();

        std::cout << "local change = " << new_action - old_action << "; total change = " << new_total - old_total << std::endl;
    }
    std::cout << "x" << std::endl;
    {
        int i = _index(gen);
        int j = _index(gen);
        int k = _index(gen);
        int l = _index(gen);

        su2 delta = su2::generate(gen);
        su2& link = lattice(i,j,k,l)->forward_x;

        double old_action, new_action;
        double old_total, new_total;

        old_action = lattice.x_link_action(i,j,k,l);
        old_total = lattice.total_action();

        link *= delta;

        new_action = lattice.x_link_action(i,j,k,l);
        new_total = lattice.total_action();

        std::cout << "local change = " << new_action - old_action << "; total change = " << new_total - old_total << std::endl;
    }
    std::cout << "y" << std::endl;
    {
        int i = _index(gen);
        int j = _index(gen);
        int k = _index(gen);
        int l = _index(gen);

        su2 delta = su2::generate(gen);
        su2& link = lattice(i,j,k,l)->forward_y;

        double old_action, new_action;
        double old_total, new_total;

        old_action = lattice.y_link_action(i,j,k,l);
        old_total = lattice.total_action();

        link *= delta;

        new_action = lattice.y_link_action(i,j,k,l);
        new_total = lattice.total_action();

        std::cout << "local change = " << new_action - old_action << "; total change = " << new_total - old_total << std::endl;
    }
    std::cout << "z" << std::endl;
    {
        int i = _index(gen);
        int j = _index(gen);
        int k = _index(gen);
        int l = _index(gen);

        su2 delta = su2::generate(gen);
        su2& link = lattice(i,j,k,l)->forward_z;

        double old_action, new_action;
        double old_total, new_total;

        old_action = lattice.z_link_action(i,j,k,l);
        old_total = lattice.total_action();

        link *= delta;

        new_action = lattice.z_link_action(i,j,k,l);
        new_total = lattice.total_action();

        std::cout << "local change = " << new_action - old_action << "; total change = " << new_total - old_total << std::endl;
    }

#undef _r
#undef _4r
}
