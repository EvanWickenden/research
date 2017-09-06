#ifndef __LATTICE_2_H__
#define __LATTICE_2_H__


#include "su2.h"
#include "callback/callback.h"

struct _Lattice
{
    struct LatticeSite
    {
        su2 links[4];   
    };

    template <int length, int width>
    struct WilsonLoop
    {
        su2 top[width], bottom[width], left[length], right[length];

        template <int mu, int nu>
        static WilsonLoop construct(const _Lattice& lattice, int t, int x, int y, int z);
        su2 contract();
    };

    virtual const _Lattice::LatticeSite *operator () (int t, int x, int y, int z) const = 0; 
    virtual _Lattice::LatticeSite *operator () (int t, int x, int y, int z) = 0; 

    virtual const int T() const = 0;
    virtual const int X() const = 0;
    virtual const int Y() const = 0;
    virtual const int Z() const = 0;

    template <int mu, int nu> su2 plaquette(const int coords[4]) const; 
    template <int mu, int nu> su2 wilson_loop(const int width, const int height, const int coords[4]) const;
    template <int mu> double link_action(const int coords[4]) const;

};


template <int T, int X, int Y, int Z>
struct Lattice : public _Lattice
{
    _Lattice::LatticeSite sites[T][X][Y][Z];

    Lattice(std::mt19937& generator, double width);
    Lattice(Callback<su2>&&);

    const _Lattice::LatticeSite *operator () (int t, int x, int y, int z) const; 
    _Lattice::LatticeSite *operator () (int t, int x, int y, int z);  

//    template <int n> static constexpr int dim(); 

    double total_action();
    void re_unitarize(); 
};


// implementation

template <int width, int length>
template <int mu, int nu>
_Lattice::WilsonLoop<width, length> _Lattice::WilsonLoop<width, length>::construct(const _Lattice& lattice, int t, int x, int y, int z)
{
    static_assert(mu != nu, "invalid wilson loop: mu == nu");

    int& _mu = (mu == 0) ? t : (mu == 1) ? x : (mu == 2) ? y : z;
    int& _nu = (nu == 0) ? t : (nu == 1) ? x : (nu == 2) ? y : z;

    WilsonLoop loop;

    for (int i = 0; i < length; i++)
    {
        loop.bottom[i] = lattice(t,x,y,z)->links[mu];
        ++_mu;
    }

    for (int i = 0; i < width; i++)
    {
        loop.right[i] = lattice(t,x,y,z)->links[nu];
        ++_nu;
    }

    for (int i = 0; i < length; i++)
    {
        --_mu;
        loop.top[i] = lattice(t,x,y,z)->links[mu].inverse();
    }

    for (int i = 0; i < width; i++)
    {
        --_nu;
        loop.left[i] = lattice(t,x,y,z)->links[nu].inverse();
    }
}

su2 _Lattice::contract()
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

#define for_each(i,j,k,l)           \
    for (int i = 0; i < T; i++)     \
    for (int j = 0; j < X; j++)     \
    for (int k = 0; k < Y; k++)     \
    for (int l = 0; l < Z; l++) 


#define __TEMPLATE__  template <int T, int X, int Y, int Z> 

__TEMPLATE__
Lattice<T,X,Y,Z>::Lattice(std::mt19937& generator, double width)
{
    for_each(i,j,k,l) 
    {
        sites[i][j][k][l].links[0] = su2::generate2(generator, width);
        sites[i][j][k][l].links[1] = su2::generate2(generator, width);
        sites[i][j][k][l].links[2] = su2::generate2(generator, width);
        sites[i][j][k][l].links[3] = su2::generate2(generator, width);
    }
}

__TEMPLATE__
Lattice<T,X,Y,Z>::Lattice(Callback<su2>&& callback)
{
    for_each(i,j,k,l)
    {
        sites[i][j][k][l].links[0] = callback();
        sites[i][j][k][l].links[1] = callback();
        sites[i][j][k][l].links[2] = callback();
        sites[i][j][k][l].links[3] = callback();
    }
}


__TEMPLATE__
const _Lattice::LatticeSite *
Lattice<T,X,Y,Z>::operator () (int t, int x, int y, int z) const
{
    return &sites[(t+T)%T][(x+X)%X][(y+Y)%Y][(z+Z)%Z];
}

__TEMPLATE__
typename _Lattice::LatticeSite *
Lattice<T,X,Y,Z>::operator () (int t, int x, int y, int z)
{
    return &sites[(t+T)%T][(x+X)%X][(y+Y)%Y][(z+Z)%Z];
}

__TEMPLATE__
void Lattice<T,X,Y,Z>::re_unitarize()
{
    for_each(i,j,k,l)
    {
        sites[i][j][k][l].links[0].re_unitarize();
        sites[i][j][k][l].links[1].re_unitarize();
        sites[i][j][k][l].links[2].re_unitarize();
        sites[i][j][k][l].links[3].re_unitarize();
    }
}


/* U_mn(n) = U_mu(n) U_nu(n + mu) U_mu^dagger(n + nu) U_n^dagger (n)  */
 
template <int mu, int nu>
su2 _Lattice::plaquette(const int coords[4]) const
{
    static_assert(0 <= mu && mu < 4 && 0 <= nu && nu < 4,  
            "invalid plaquette arguments: out of bounds");
    static_assert(mu != nu, "invalid plaquette arguments: mu == nu");

    static int *copy;
    copy = const_cast<int *>(coords);

    const LatticeSite& bottom_left = (*this)(copy);

    ++copy[mu];
    const LatticeSite& bottom_right = (*this)(copy);
    --copy[mu]; 

    ++copy[nu];
    const LatticeSite& top_left = (*this)(copy);
    --copy[nu];

    return bottom_left.links[mu] * bottom_right.links[nu]
        * top_left.links[mu].inverse() * bottom_left.links[nu].inverse();
}

template <int mu, int nu>
su2 _Lattice::wilson_loop(int width, int height, const int coords[4]) const
{
    static su2 product; 
    static int i;
    static int *copy;

    product = su2::identity();
    copy = const_cast<int *>(coords); 

    for (i = 0; i < width; i++)
    {
        product *= (*this)(copy).links[mu];
        ++copy[mu];
    }
    for (i = 0; i < height; i++)
    {
        product *= (*this)(copy).links[nu];
        ++copy[nu];
    }
    for (i = 0; i < width; i++)
    {
        --copy[mu];
        product *= (*this)(copy).links[mu].inverse();
    }
    for (i = 0; i < height; i++)
    {
        --copy[nu];
        product *= (*this)(copy).links[nu].inverse();
    }
    return product;
}




__TEMPLATE__
template <int mu>
double Lattice<T,X,Y,Z>::link_action(const int coords[4]) const
{
    static_assert(0 <= mu && mu < 4, "invalid link_action arguments: out of bounds");

    constexpr int a = (mu+1)%4; 
    constexpr int b = (mu+2)%4; 
    constexpr int c = (mu+3)%4;

    constexpr int offset_a = Lattice<T,X,Y,Z>::dim<a>() - 1;
    constexpr int offset_b = Lattice<T,X,Y,Z>::dim<b>() - 1;
    constexpr int offset_c = Lattice<T,X,Y,Z>::dim<c>() - 1;

    static double sum;
    static int *copy;

    copy = const_cast<int *>(coords);

    sum = plaquette<mu, a>(coords).half_trace()
        + plaquette<mu, b>(coords).half_trace()
        + plaquette<mu, c>(coords).half_trace();

    copy[a] += offset_a; 
    sum += plaquette<mu, a>(coords).half_trace();
    copy[a] -= offset_a;

    copy[b] += offset_b;
    sum += plaquette<mu, b>(coords).half_trace();
    copy[b] -= offset_b;

    copy[c] += offset_c;
    sum += plaquette<mu, c>(coords).half_trace();
    copy[c] -= offset_c;

    return -sum;
}

__TEMPLATE__
double Lattice<T,X,Y,Z>::total_action()
{
    double sum = 0;
    for_each(i,j,k,l)
    {
        int coords[] = {i,j,k,l};
        sum += plaquette<0,1>(coords).half_trace()
        + plaquette<0,2>(coords).half_trace()
        + plaquette<0,3>(coords).half_trace()
        + plaquette<1,2>(coords).half_trace()
        + plaquette<1,3>(coords).half_trace()
        + plaquette<2,3>(coords).half_trace();
    }
    return -sum;
}

#undef __TEMPLATE__
#undef for_each

#endif
