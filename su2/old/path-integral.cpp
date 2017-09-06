
#include "path-integral.h"

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
    su2& operator[](int i)
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


template <int size>
static inline void 
update_configuration(Lattice& lattice, double beta, 
        std::mt19937& generator, std::uniform_real_distribution<double>& _accept, 
        std::uniform_int_distribution<int>& _matrix_index, RandomMatrices<size>& random_matrices)
{
    int i, j, k, l;
    double old_action, new_action;
    su2 copy, delta;

    auto propose_change = [&] (su2& link, double (Lattice::*action)(int,int,int,int))
    {
        delta = random_matrices[_matrix_index(generator)];

        old_action = (lattice.*action)(i,j,k,l);
        copy = link;
        link = delta * link;
        new_action = (lattice.*action)(i,j,k,l);

        lattice.acceptance_ratio[1]++;

        if (new_action > old_action && 
                _accept(generator) > exp( 0.5 * beta * (old_action - new_action))) 
            link = copy;
        else
            lattice.acceptance_ratio[0]++;
    };

    for (i = 0; i < lattice.dimensions[0]; i++)
    {
     for (j = 0; j < lattice.dimensions[1]; j++)
     {
      for (k = 0; k < lattice.dimensions[2]; k++)
      {
       for (l = 0; l < lattice.dimensions[3]; l++)
       {
           Lattice::LatticeSite *site = lattice(i,j,k,l);
           propose_change(site->forward_t, &Lattice::t_link_action);
           propose_change(site->forward_x, &Lattice::x_link_action);
           propose_change(site->forward_y, &Lattice::y_link_action);
           propose_change(site->forward_z, &Lattice::z_link_action);
       }
      }
     }
    }
}


void path_integral(Lattice& lattice, double beta, int nr_configurations, double width, Observable& observable)
{
    std::random_device ran;
    std::mt19937 generator(ran());

    std::uniform_real_distribution<double> _accept(0, 1);
    std::uniform_int_distribution<int> _matrix_index(0, 63);
    RandomMatrices<64> random_matrices(generator, width);   

    int k = 0;

    for (int i = 0; i < nr_configurations; i++)
    {
        update_configuration<64>(lattice, beta, generator, _accept, _matrix_index, random_matrices);
        observable(lattice);

        if (k++ % 10 == 0)
            lattice.re_unitarize();
    }
}


