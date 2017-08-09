
#include <iostream>

#include "callback.h"


struct Test : public Callback<double, double>
{
    double operator () (double x)
    {
        return 4;
    }
};


typedef Callback<void, double> O;

