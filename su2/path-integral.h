#ifndef __PATH_INTEGRAL_H__
#define __PATH_INTEGRAL_H__

#include <random>

#include "lattice.h"
#include "callback/callback.h"
#include "su2.h"

typedef Callback<void, Lattice&> Observable;

void path_integral(Lattice& lattice, double beta, int nr_configurations, double width, Observable& o);





#endif
