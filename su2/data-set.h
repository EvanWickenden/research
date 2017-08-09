
#include <vector>
#include <iostream>
#include <cstdlib>
#include <math.h>

struct DataSet
{
    double *data;
    int n;
    int offset;
    int index;
    double mean, integrated_autocorrelation_time, exponential_autocorrelation_time, std_dev;
    double c0;

    DataSet(int n);
    ~DataSet();

    DataSet& operator = (const DataSet& rhs);

    void shift(int difference);
    void store(double datum);

    void _mean();
    double _autocorrelation_fn(int t);
    void _c0();
    void _exponential_autocorrelation_time(double smoothing = 0);
    void _integrated_autocorrelation_time(int window = 1, double a = 8);
    void _std_dev();

    void analyze();

    friend std::ostream& operator << (std::ostream& stream, DataSet& data);
};
