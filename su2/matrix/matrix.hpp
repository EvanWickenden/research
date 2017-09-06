
#ifndef __MATRIX_H__
#define __MATRIX_H__


#include "tuple.h"


template <int N>
struct Matrix
{
#define for_each(_index)                        \
    for (int _index = 0; _index < N; _index++)

    Complex components[N][N];

    Complex& operator () (int i, int j) { return components[i][j]; }
    const Complex& operator () (int i, int j) const { return components[i][j]; }

    Tuple<N> row(int i) const
    { 
        Tuple<N> t; 
        for_each(j) t[j] = components[i][j];
        return t;
    }
    Tuple<N> column(int i) const
    {
        Tuple<N> t; 
        for_each(j) t[j] = components[j][i];
        return t;
    }

    friend std::ostream& operator << (std::ostream& stream, const Matrix& M)
    {
        for_each(i) stream << "|  " << M.row(i) << "  |" << std::endl;
        return stream;
    }

    static Matrix<N> identity()
    {
        Matrix<N> M;
        for_each(i) M.components[i][i] = 1;
        return M;
    }


    Matrix& operator *= (const double c)
    {
        for_each(i) for_each(j) components[i][j] *= c;
        return *this;
    }
    Matrix& operator *= (const Matrix& M)
    {
        Matrix copy = *this;
        for_each(i) for_each(j) components[i][j] = copy.row(i) * M.column(j);
        return *this;
    }
    Matrix& operator += (const Matrix& M)
    {
        for_each(i) for_each(j) components[i][j] += M(i,j);
        return *this;
    }
};


template <int N>
int operator == (const Matrix<N>& lhs, const Matrix<N>& rhs)
{
    for_each(i) for_each(j)
        if (lhs(i,j) != rhs(i,j)) return 0;
    return 1;
}

template <int N>
Matrix<N> operator * (const double c, const Matrix<N>& rhs)
{
    return Matrix<N>(rhs) *= c;
}
template <int N>
Matrix<N> operator * (const Matrix<N>& lhs, const Matrix<N>& rhs)
{
    return Matrix<N>(lhs) *= rhs;
}
template <int N>
Matrix<N> operator + (const Matrix<N>& lhs, const Matrix<N>& rhs)
{
    return Matrix<N>(lhs) += rhs;
}

template <int N>
Matrix<N> exp(const Matrix<N>& M, int precision = 10)
{
    double c = 1;
    Matrix<N> product = Matrix<N>::identity();
    Matrix<N> sum = Matrix<N>::identity();

    for (int i = 1; i < precision; i++)
    {
        product *= M;
        c /= i;
        sum +=  c * product;
    }

    return sum;
}

#undef for_each

#endif 

