
#include <random>
#include <iostream>

#include "path-integral.h"
#include "lattice.h"
#include "data-set.h"


struct _Quark 
{
    int R, T;
    DataSet data;

    _Quark(int R, int T, int n) :
        R(R),
        T(T),
        data(n)
    { }

    void operator () (Lattice& lattice)
    {
        su2 product = su2::identity();
        int x = 0, t = 0;
        for (; t < T; t++)
            product *= lattice(t,x,0,0)->forward_t.inverse();
        for (; x < R; x++)
            product *= lattice(t,x,0,0)->forward_x.inverse();
        for (; t > 0; t--)
            product *= lattice(t,x,0,0)->forward_t;
        for (; x > 0; x--)
            product *= lattice(t,x,0,0)->forward_x;
        data.store(product.trace().real());
    }

    void analyze()
    {
        data.analyze();
    }

    friend std::ostream& operator << (std::ostream& stream, _Quark& q)
    {
        stream << "quark path - T = " << q.T << ", R = " << q.R << ";\n";
        stream << q.data;
        return stream;
    }
};

/*
template <int R, int T>
struct Quarks : public Observable
{
    _Quark* quarks[R][T];  

#define for_each(i,j)                   \
    for (int i = 0; i < R; i++)         \
    for (int j = 0; j < T; j++)

    Quarks(int n) 
    {
        for_each(i,j) 
            quarks[i][j] = new _Quark(i+1, j+1, n);
    }

    ~Quarks() 
    {
        for_each(i,j) 
            delete quarks[i][j];
    }

    void operator () (Lattice& lattice)
    {
        for_each(i,j) 
            (*quarks[i][j])(lattice);
    }

    void analyze()
    {
        for_each(i,j) 
            quarks[i][j]->analyze();
    }

    friend std::ostream& operator << (std::ostream& stream, Quarks& q)
    {
        for_each(i,j) 
            stream << *q.quarks[i][j];
        return stream;
    }
#undef for_each
};
*/

template <int R, int T>
struct Quarks : public Observable
{
    _Quark* quarks[R];

#define for_each(i)                 \
    for (int i = 0; i < R; i++)

    Quarks(int n)
    {
        for_each(i) quarks[i] = new _Quark(i+1, i+1, n);
    }

    ~Quarks()
    {
        for_each(i) delete quarks[i];
    }

    void operator () (Lattice& lattice)
    {
        for_each(i) (*quarks[i])(lattice);
    }

    void analyze()
    {
        for_each(i) quarks[i]->analyze();
    }

    friend std::ostream& operator << (std::ostream& stream, Quarks& q)
    {
        for_each(i) stream << *q.quarks[i] << std::endl;
        return stream;
    }
#undef for_each
};

struct Nil : public Observable
{
    void operator () (Lattice& lattice)
    { }
};


int main()
{
    std::random_device ran;
    std::mt19937 generator(ran());

    int lattice_dimensions[4] = {8, 8, 4, 4};
    double beta = 5;
    double width = 0.1;
    int nr_configurations = 100;

    Lattice lattice(lattice_dimensions, generator, 0.1);
    Quarks<6,6> quarks(nr_configurations);
    path_integral2(lattice, beta, nr_configurations, width, quarks, generator);

    std::cout << lattice.acceptance_ratio[0] / (double) lattice.acceptance_ratio[1] << std::endl;
    quarks.analyze();
    std::cout << quarks << std::endl;
}




