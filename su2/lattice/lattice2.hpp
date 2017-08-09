#ifndef __LATTICE_2_H__
#define __LATTICE_2_H__


#include "../su2.h"


template <int T, int X, int Y, int Z>
struct Lattice2 
{
    struct LatticeSite
    {
        su2 links[4];   
    } 
    sites[T][X][Y][Z];

    Lattice2(std::mt19937& generator, double width);

    const Lattice2::LatticeSite *operator () (int t, int x, int y, int z) const; 
    Lattice2::LatticeSite *operator () (int t, int x, int y, int z); 

    void re_unitarize(); 

    template <int mu, int nu> su2 plaquette(int coords[4]); 
    template <int n> constexpr int dim() const; 
    template <int mu> double link_action(int coords[4]);

    double total_action();
};

// implementation

#define for_each(i,j,k,l)           \
    for (int i = 0; i < T; i++) {   \
    for (int j = 0; j < X; j++) {   \
    for (int k = 0; k < Y; k++) {   \
    for (int l = 0; l < Z; l++) {

#define __end__ }}}}

template <int T, int X, int Y, int Z>
Lattice2<T,X,Y,Z>::Lattice2(std::mt19937& generator, double width)
{
    for_each(i,j,k,l)
        sites[i][j][k][l].links[0] = su2::generate2(generator, width);
        sites[i][j][k][l].links[1] = su2::generate2(generator, width);
        sites[i][j][k][l].links[2] = su2::generate2(generator, width);
        sites[i][j][k][l].links[3] = su2::generate2(generator, width);
    __end__
}


template <int T, int X, int Y, int Z>
const typename Lattice2<T,X,Y,Z>::LatticeSite *
Lattice2<T,X,Y,Z>::operator () (int t, int x, int y, int z) const
{
    return &sites[t%T][x%X][y%Y][z%Z];
}

template <int T, int X, int Y, int Z>
typename Lattice2<T,X,Y,Z>::LatticeSite *
Lattice2<T,X,Y,Z>::operator () (int t, int x, int y, int z)
{
    return &sites[t%T][x%X][y%Y][z%Z];
}

template <int T, int X, int Y, int Z>
void Lattice2<T,X,Y,Z>::re_unitarize()
{
    for_each(i,j,k,l)
        sites[i][j][k][l].links[0].re_unitarize();
        sites[i][j][k][l].links[1].re_unitarize();
        sites[i][j][k][l].links[2].re_unitarize();
        sites[i][j][k][l].links[3].re_unitarize();
    __end__
}

// U_mn(n) = U_mu(n) U_nu(n + mu) U_mu^dagger(n + nu) U_n^dagger (n) 
template <int T, int X, int Y, int Z>
template <int mu, int nu>
su2 Lattice2<T,X,Y,Z>::plaquette(int coords[4])
{
    static_assert(0 <= mu && mu < 4 && 0 <= nu && nu < 4,  
            "invalid plaquette arguments: out of bounds");
    static_assert(mu != nu, "invalid plaquette arguments: mu == nu");

#define _site(ar) \
    (sites[ar[0] % T][ar[1] % X][ar[2] % Y][ar[3] % Z])

    if constexpr (mu < nu)
    {
        LatticeSite& bottom_left = _site(coords);
        ++coords[mu];
        LatticeSite& bottom_right = _site(coords);
        --coords[mu]; ++coords[nu];
        LatticeSite& top_left = _site(coords);
        --coords[nu];

        return bottom_left.links[mu] * bottom_right.links[nu]
            * top_left.links[mu].inverse() * bottom_left.links[nu].inverse();
    }
    else
    {
        LatticeSite& bottom_left = _site(coords);
        ++coords[nu];
        LatticeSite& bottom_right = _site(coords);
        --coords[nu]; ++coords[mu];
        LatticeSite& top_left = _site(coords);
        --coords[mu];

        return bottom_left.links[nu] * bottom_right.links[mu]
            * top_left.links[nu].inverse() * bottom_left.links[mu].inverse();
    }
#undef _site
}

template <int T, int X, int Y, int Z>
template <int n>
constexpr int Lattice2<T,X,Y,Z>::dim() const
{
    if constexpr (n == 0) return T;
    else if constexpr (n == 1) return X;
    else if constexpr (n == 3) return X;
    else return Z;
}

template <int T, int X, int Y, int Z>
template <int mu>
double Lattice2<T,X,Y,Z>::link_action(int coords[4])
{
    static_assert(0 <= mu && mu < 4, "invalid link_action arguments: out of bounds");
    double sum = 0; 

    constexpr int a = (mu+1)%4; 
    constexpr int b = (mu+2)%4; 
    constexpr int c = (mu+3)%4;

    sum += plaquette<mu, a>(coords).trace().real()
        + plaquette<mu, b>(coords).trace().real()
        + plaquette<mu, c>(coords).trace().real();

    coords[a] += dim<a>() - 1;
    sum += plaquette<mu, a>(coords).trace().real();
    ++coords[a]; 

    coords[b] += dim<b>() - 1;
    sum += plaquette<mu, b>(coords).trace().real();
    ++coords[b];

    coords[c] += dim<c>() -1;
    sum += plaquette<mu, c>(coords).trace().real();
    ++coords[c];

    return sum;
}

template <int T, int X, int Y, int Z>
double Lattice2<T,X,Y,Z>::total_action()
{
    double sum = 0;
    for_each(i,j,k,l)
        int coords[] = {i,j,k,l};
        sum += plaquette<0,1>(coords).trace().real()
        + plaquette<0,2>(coords).trace().real()
        + plaquette<0,3>(coords).trace().real()
        + plaquette<1,2>(coords).trace().real()
        + plaquette<1,3>(coords).trace().real()
        + plaquette<2,3>(coords).trace().real();
    __end__
    return sum;
}

#undef for_each

#endif
