
#ifndef __TUPLE_H__
#define __TUPLE_H__

#include <complex>
#include <iostream>

typedef std::complex<double> Complex;

template <int N>
struct Tuple
{
#define for_each(_index)                        \
    int _index;                                 \
    for (_index = 0; _index < N; _index++)

    Complex components[N];

    Complex& operator[](int i) { return components[i]; }
    const Complex& operator[](int i) const { return components[i]; }

//    Tuple(){};
//    Tuple(Complex components[N]) { for_each(i) this->components[i] = components[i]; }
//    Tuple(const Tuple& t) { for_each(i) components[i] = t[i]; }
//    Tuple(const Tuple&& t) { for_each(i) components[i] = t[i]; }
//    Tuple& operator = (const Tuple& t) { for_each(i) components[i] = t[i]; return *this; }
//    Tuple& operator = (const Tuple&& t) { for_each(i) components[i] = t[i]; return *this; }

    friend std::ostream& operator << (std::ostream& stream, Tuple<N> t)
    {
        stream << t[0];
        for (int i = 1; i < N; i++)
            stream << ", " << t[i];
        return stream;
    }

    /* dot product */
    Complex operator *= (const Tuple& rhs)
    {
        Complex sum = 0;
        for_each(i) sum += components[i] * rhs[i];
        return sum;
    }

    Tuple& operator *= (const Complex c)
    {
        for_each(i) components[i] *= c;
        return *this;
    }

    Tuple& operator += (const Tuple& rhs)
    {
        for_each(i) components[i] += rhs[i];
        return *this;
    }
};

template <int N>
Complex operator * (const Tuple<N>& lhs, const Tuple<N>& rhs)
{
    return Tuple<N>(lhs) *= rhs;
}


#undef for_each



#endif 
