
#include <math.h>

#include "path-integral.h"


struct X1X2 : public Observable
{
	int i, j;
	double data;

	X1X2(int i, int j) : i(i), j(j), data(0) { }

	void operator()(const Lattice& lattice)
	{
		data += lattice[i] * lattice[j];
	}
};

struct Correlation : public Observable
{
	int index, offset;
	int count;
	double *data;

	Correlation(int index, int offset) :
		index(index), offset(offset), count(0),
		data(new double[offset])
	{
		int i;
		for (i = 0; i < offset; i++)
			data[i] = 0;
	}
	~Correlation(){ delete[] data; }

	void operator()(const Lattice& lattice)
	{
		int i;
		for (i = 0; i < offset; i++)
		{
			data[i] += lattice[index] * lattice[index + i];
		}
		count++;
	}

	double operator[](int i) { return data[i] / (float) count; }

	friend std::ostream& operator<<(std::ostream& s, Correlation& c)
	{
		int i;
		s << "correlation measurement starting at index " << c.index << std::endl;
		for (i = 0; i < c.offset; i++)
		{
			s << c[i] << "\n";
		}
		return s;
	}
};


struct HO : public DeltaAction
{
	double tau;
	double alpha;

	HO(double tau, double omega) : tau(tau), alpha(0.5*tau*tau*omega*omega) {}

	double operator()(const Lattice& lattice, double Delta_x, int index) 
	{
		return (2 / tau) * Delta_x * ((1 + alpha) * (2 * lattice[index] + Delta_x) - lattice[index + 1] - lattice[index - 1]);
	}
};


int main()
{
	const long steps = 500000000;
	const double tau = 1 / (float) 40;
	const double omega = 3;
	const int lattice_size = 200;

	/* analytic solution for 2 elt lattice */
//	double prediction = 1 / ( tau * omega * omega * (4 + tau * tau * omega * omega ) );
//	std::cout << "prediction = " << prediction << std::endl;

	Lattice lattice(lattice_size);
	Correlation cor(20, lattice_size + 1);
	HO ho(tau, omega);

	int i;
	for (i = 0; i < lattice_size; i++)
	{
		lattice[i] = 0;
	}

	path_integral(lattice, steps, 0.09, ho, cor);

	std::cout << cor;
}
