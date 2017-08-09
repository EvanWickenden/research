
#include <random>
#include <iostream>
#include <assert.h>

#include "../lattice.h"
#include "../su2.h"

#define _id su2(1,0)
#define SUCCESS()   std::cout << "success!" << std::endl;

#define _r  su2::generate(gen)
#define _4r _r, _r, _r, _r


int main()
{
    std::random_device ran;
    std::mt19937 gen(ran());

    {
        Lattice lattice((int[]) {2, 2, 1, 1}, gen);
        su2 plaquette = lattice.tx_plaquette(0,0,0,0);
        double det = std::norm(plaquette.a) + std::norm(plaquette.b);
        std::cout << "randomly generated plaquette determinant = " << det << std::endl;
        SUCCESS();
    }

    {
        Lattice lattice((int[]) {1,1,1,1}, gen);

        std::cout << lattice.tx_plaquette(0,0,0,0) << std::endl;

    }

    Lattice::LatticeSite sites[2][2] =
    {
        {
            { _4r },
            { _4r }
        },
        {
            { _4r },
            { _4r }
        }
    };

    {
        su2 plaquette = sites[0][0].forward_t.inverse()
            * sites[1][0].forward_x.inverse()
            * sites[0][1].forward_t
            * sites[0][0].forward_x;

        Lattice lattice( (int[]) {2,2,1,1}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.tx_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: tx" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }
    {
        su2 plaquette = sites[0][0].forward_t.inverse()
            * sites[1][0].forward_y.inverse()
            * sites[0][1].forward_t
            * sites[0][0].forward_y;

        Lattice lattice( (int[]) {2,1,2,1}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.ty_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: ty" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }
    {
        su2 plaquette = sites[0][0].forward_t.inverse()
            * sites[1][0].forward_z.inverse()
            * sites[0][1].forward_t
            * sites[0][0].forward_z;

        Lattice lattice( (int[]) {2,1,1,2}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.tz_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: tz" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }
    {
        su2 plaquette = sites[0][0].forward_x.inverse()
            * sites[1][0].forward_y.inverse()
            * sites[0][1].forward_x
            * sites[0][0].forward_y;

        Lattice lattice( (int[]) {1,2,2,1}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.xy_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: xy" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }
    {
        su2 plaquette = sites[0][0].forward_x.inverse()
            * sites[1][0].forward_z.inverse()
            * sites[0][1].forward_x
            * sites[0][0].forward_z;

        Lattice lattice( (int[]) {1,2,1,2}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.xz_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: xz" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }
    {
        su2 plaquette = sites[0][0].forward_y.inverse()
            * sites[1][0].forward_z.inverse()
            * sites[0][1].forward_y
            * sites[0][0].forward_z;

        Lattice lattice( (int[]) {1,1,2,2}, (Lattice::LatticeSite *) sites);
        su2 p2 = lattice.yz_plaquette(0,0,0,0);

        std::cout << "testing plaquette link product: yz" << std::endl;
        assert(plaquette == p2);
        SUCCESS();
    }


    Lattice::LatticeSite sites2[2][2][2][2] = 
    {
        {
            {
                {
                    { _4r },
                    { _4r }
                },
                {
                    { _4r },
                    { _4r }
                }
            }
        }
    };

}
