
#include "callback/callback.h"
#include "lattice/lattice2.hpp"
#include "path-integral2.hpp"

#include "data-set.h"

#define abs_val(a) ({ double _a = a; (_a < 0) ? -_a : _a; })

static const int T = 8, X = 8, Y = 8, Z = 8;

struct _Quark : public Observable<T,X,Y,Z>
{
    int length, time;
    DataSet data;

    _Quark(int length, int time, int n) :
        length(length),
        time(time),
        data(n)
    { }

    void operator () (const Lattice2<T,X,Y,Z>& lattice) override
    {
        su2 product = su2::identity();
        int t = 0, x = 0, y = 0, z = 0;
        for (; t < time; t++)
            product *= lattice(t,x,y,z)->links[0];
        for (; x < length; x++)
            product *= lattice(t,x,y,z)->links[1];
        for (; t > 0; t--)
            product *= lattice(t,x,y,z)->links[0].inverse();
        for (; x > 0; x--)
            product *= lattice(t,x,y,z)->links[1].inverse();
        data.store(product.trace().real());
    }

    void analyze()
    {
        data.analyze();
    }

    friend std::ostream& operator << (std::ostream& stream, _Quark& q)
    {
        stream << "quark path - T = " << q.time << ", R = " << q.length << ";\n";
        stream << q.data;
        return stream;
    }
};

template <int length>
struct Quarks : public Observable<T,X,Y,Z>
{
    _Quark* quarks[length];

#define for_each(i)                 \
    for (int i = 0; i < length; i++)

    Quarks(int n)
    {
        for_each(i) quarks[i] = new _Quark(i+1, 6, n);
    }

    ~Quarks()
    {
        for_each(i) delete quarks[i];
    }

    void operator () (const Lattice2<T,X,Y,Z>& lattice) override
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

int main()
{
    double beta = 4;
    double nr_configurations = 100000;
    double width = 0.15;

    Quarks<7> q(nr_configurations);

    std::random_device ran;
    std::mt19937 gen(ran());

    static Lattice2<T,X,Y,Z> lattice(gen, 0.05);
    
    path_integral2<T,X,Y,Z>(lattice, beta, nr_configurations, width, q, gen);

    q.analyze();
    std::cout << q << std::endl;
}
