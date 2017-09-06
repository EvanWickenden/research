
#include <su2.h>
#include <lattice.hpp>

template <int length, int width>
struct WilsonLoop
{
    su2 top[width], bottom[width], left[length], right[length];

    template <int mu, int nu>
    static WilsonLoop construct(const _Lattice& lattice, int t, int x, int y, int z)
    {
        static_assert(mu != nu, "invalid wilson loop: mu == nu");

        int& _mu = (mu == 0) ? t : (mu == 1) ? x : (mu == 2) ? y : z;
        int& _nu = (nu == 0) ? t : (nu == 1) ? x : (nu == 2) ? y : z;

        WilsonLoop loop;
        _Lattice::LatticeSite *site;

        for (int i = 0; i < length; i++)
        {
            bottom[l] = lattice(t,x,y,z).links[mu];
            ++_mu;
        }

        for (int i = 0; i < width; i++)
        {
            right[i] = lattice(t,x,y,z).links[nu];
            ++_nu;
        }

        for (int i = 0; i < length; i++)
        {
            --_mu;
            top[i] = lattice(t,x,y,z).links[mu].inverse();
        }

        for (int i = 0; i < width; i++)
        {
            --_nu;
            left[i] = lattice(t,x,y,z).links[nu].inverse();
        }
    }

    su2 contract()
    {
        su2 product = su2::identity();

        for (int i = 0; i < length; i++)
            product *= bottom[i];
        for (int i = 0; i < width; i++)
            product *= right[i];
        for (int i = 0; i < length; i++)
            product *= top[i];
        for (int i = 0; i < width; i++)
            product *= left[i];

        return product;
    }
};

