
#include <vector>
#include <iostream>
#include <cstdlib>
#include <math.h>

#include "data-set.h"

#define min(a, b) ({ double _a = a; double _b = b; (_a < _b) ? _a : _b; })


DataSet::DataSet(int n) :
    data(NULL),
    n(n),
    offset(0),
    index(0)
{
    data = (double *) malloc(sizeof(double) * n);
    if (data == NULL)
    {
        std::cout << "malloc() failed" << std::endl;
        exit(1);
    }
}

DataSet::~DataSet() 
{ 
    free(data - offset); 
}


void DataSet::shift(int difference)
{
    offset += difference;
    data += difference;
    n -= difference;
}

void DataSet::store(double datum)
{
    data[index++] = datum;
}

void DataSet::_mean()
{
    double sum = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        sum += data[i];
    }
    mean = sum / n;
}


double DataSet::_autocorrelation_fn(int t)
{
    /* TODO: remove error checking */
    if (t >= n)
    {
        std::cout << "autocorrelation_fn out of bounds" << std::endl;
        exit(1);
    }

    double sum = 0;
    int i; 
    for (i = 0; i < n - t; i++)
    {
        sum += (data[i] - mean) * (data[i + t] - mean);
    }

    return sum / (n - t); 
}

void DataSet::_c0()
{
    c0 = _autocorrelation_fn(0);
}

void DataSet::_exponential_autocorrelation_time(double smoothing)
{
    double c = 0;
    double value = 0;
    int index = n / 2;
    double max = 0;
    double min = 0;

    c = c0;

#define _next_val                                                                   \
    if (index >= n) goto done;                                                      \
    c = (1 - smoothing) * _autocorrelation_fn(index) + smoothing * c;               \
    value = -index / log(c / c0);                                                   \
    index++

increasing:
    _next_val;
    if (value < max)
    {
        min = value;
        goto decreasing;
    }
    max = value;
    goto increasing;

decreasing:
    _next_val;
    if (value > min)
    {
        max = value;
        goto increasing;
    }
    min = value;
    goto decreasing;

done:
    exponential_autocorrelation_time = max;
#undef _next_val
}

void DataSet::_integrated_autocorrelation_time(int window, double a)
{
    double tau_int = 0.5;
    int i = 1;

    while (i < n
            && i < window)
    {
        tau_int += fabs(_autocorrelation_fn(i++) / c0);
    }
    while (i < n
            && window < a*tau_int)
    {
        tau_int += fabs(_autocorrelation_fn(i++) / c0);
        window++;
    }
    integrated_autocorrelation_time = tau_int;
}

void DataSet::_std_dev()
{
    std_dev = sqrt(2 * integrated_autocorrelation_time * c0 / n);
}

void DataSet::analyze()
{
    shift(0.2 * n);
    _mean();
    _c0();
//    _exponential_autocorrelation_time(0.3);
//    shift(min(exponential_autocorrelation_time, 0.2 * n));
//    _mean();
//    _c0();
    _integrated_autocorrelation_time();
    _std_dev();
}

std::ostream& operator<<(std::ostream& stream, DataSet& data)
{
    stream << "mean = " << data.mean << ", exponential autocorrelation time = " << data.exponential_autocorrelation_time
        << ", integrated autocorrelation time = " << data.integrated_autocorrelation_time 
        << ", standard deviation = " << data.std_dev << std::endl;
    return stream;
}
