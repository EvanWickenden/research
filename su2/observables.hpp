#ifndef __OBSERVABLES_HPP__
#define __OBSERVABLES_HPP__

#include "data-set.h"
#include "lattice.hpp"
#include "su2.h"
#include "path-integral.hpp"

template <int T, int X, int Y, int Z>
struct AvPlaquette : public Observable<T,X,Y,Z>, public Analyzable
{
    DataSet d[6];

    AvPlaquette(int n)
        : d{n,n,n,n,n,n}
    {}

    void operator () (const Lattice<T,X,Y,Z>& l)
    {
        int coords[4] = {1,1,1,1};
        d[0].store( l.template wilson_loop<0,1>(2,2,coords).half_trace() );
        d[1].store( l.template wilson_loop<0,2>(2,2,coords).half_trace() );
        d[2].store( l.template wilson_loop<0,3>(2,2,coords).half_trace() );
        d[3].store( l.template wilson_loop<1,2>(2,2,coords).half_trace() );
        d[4].store( l.template wilson_loop<1,3>(2,2,coords).half_trace() );
        d[5].store( l.template wilson_loop<2,3>(2,2,coords).half_trace() );
    }


#define for_each(i)    for (int i = 0; i < 6; i++)

    void analyze()
    {
        for_each(i) d[i].analyze();
    }

    void write(std::ostream& stream)
    {
        for_each(i) stream << d[i] << std::endl;
    }

#undef for_each
};


#endif
