#ifndef __PATH_INTEGRAL_H__
#define __PATH_INTEGRAL_H__

#include "warray.h"

typedef WArray<double> Lattice;

struct Observable
{
	virtual void operator()(const Lattice& lattice) = 0;
};

struct DeltaAction
{
	virtual double operator()(const Lattice& lattice, double Delta_x, int index) = 0;
};

void path_integral(Lattice& lattice, long steps, double delta_x_width, DeltaAction& DS, Observable& observable);

#endif
