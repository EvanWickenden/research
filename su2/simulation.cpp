
#include <assert.h>

#include "callback/callback.h"
#include "lattice.hpp"
#include "path-integral.hpp"
#include "data-set.h"
#include "observables.hpp"

#define abs_val(a) ({ double _a = a; (_a < 0) ? -_a : _a; })

static const int T = 3, X = 8, Y = 8, Z = 8;

struct _Quark 
{
    int width, height;
    DataSet data;

    _Quark(int width, int height, int n) :
        width(width),
        height(height),
        data(n)
    { }

    void operator () (const Lattice<T,X,Y,Z>& lattice)
    {
        data.store(lattice.wilson_loop<0,1>(width, height, (int[4]){0,0,0,0}).half_trace());
    }

    void analyze()
    {
        data.analyze();
    }

    friend std::ostream& operator << (std::ostream& stream, _Quark& q)
    {
        stream << "quark path - T = " << q.height << ", R = " << q.width << ";\n";
        stream << q.data;
        return stream;
    }
};

template <int length>
struct Quarks : public Observable<T,X,Y,Z>
{
    _Quark* quarks[length][length];

    struct Data
    {
        double mean, std_dev;
        int T, R;
    };

    std::vector<Data> data[(length + 1) * (length + 1)];

#define for_each(i,j)                       \
    for (int i = 0; i < length; i++)        \
    for (int j = 0; j < length; j++)

    Quarks(int n) 
    {
        for_each(i,j) quarks[i][j] = new _Quark(i+1, j+1, n);
    }

    ~Quarks()
    {
        for_each(i,j) delete quarks[i][j];
    }

    void operator () (const Lattice<T,X,Y,Z>& lattice) override
    {
        for_each(i,j) (*quarks[i][j])(lattice);
    }

    void analyze()
    {
        for_each(i,j) 
        {
            quarks[i][j]->analyze();
            data[(i+1)*(j+1)].push_back((Data){quarks[i][j]->data.mean, quarks[i][j]->data.std_dev, i, j});
        }
    }

    friend std::ostream& operator << (std::ostream& stream, Quarks& q)
    {
        for (int i = 1; i <= length * length; i++)
        {
            std::cout << "Area: " << i << std::endl;
            for (auto _data = q.data[i].begin(); _data != q.data[i].end(); ++_data)
            {
                std::cout << "[" << _data->mean << ", " << _data->std_dev << ", " << _data->T << ", " << _data->R << "] ";
            }
            std::cout << "\n" << std::endl;
        }

        return stream;
    }

#undef for_each
};

struct Nil : public Observable<T,X,Y,Z>
{
    void operator () (const Lattice<T,X,Y,Z>& lattice) override {}
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
    Nil n;

    path_integral<T,X,Y,Z>(lattice, beta, nr_configurations, width, n, gen);

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

