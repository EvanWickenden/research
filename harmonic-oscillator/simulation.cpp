
#include <math.h>

#include "path-integral.h"
#include "data.h"

/* new functions */

struct Measurement
{
    Observable &o;
    Measurement(Observable& o) : o(o) {};
};


/*    */



struct Correlation : public Observable
{
	int index, offset;
	int count;
	Data **data_sets;

	Correlation(int index, int offset, int data_set_size) :
		index(index), offset(offset), count(0),
        data_sets(new Data *[offset])
	{
		int i;
		for (i = 0; i < offset; i++)
        {
           data_sets[i] = new Data(data_set_size);
        }
	}
	~Correlation()
    { 
        int i;
        for (i = 0; i < offset; i++)
        {
            delete data_sets[i];
        }
        delete[] data_sets;
    }

	void operator()(const Lattice& lattice)
	{
		int i;
		for (i = 0; i < offset; i++)
			data_sets[i]->data[count] = lattice[index] * lattice[index + i];
        count++;
	}

    void analyze_data_sets()
    {
        int i;
        for (i = 0; i < offset; i++)
            data_sets[i]->analyze_data_set();
    }

	friend std::ostream& operator<<(std::ostream& s, Correlation& c)
	{
		int i;
        s << "\n means: " << std::endl;
		for (i = 0; i < c.offset; i++)
		{
            s << c.data_sets[i]->mean << std::endl;
		}
        s << "\n std devs: " << std::endl;
		for (i = 0; i < c.offset; i++)
		{
            s << c.data_sets[i]->std_dev << std::endl;
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
	const double tau = 0.2;
	const double omega = 1;
	const int lattice_size = 100;

//    std::cout << 1 / (tau * omega * omega * 4) << std::endl;

	Lattice lattice(lattice_size);
	Correlation cor(0, lattice_size, steps / lattice_size);
	HO ho(tau, omega);

	int i;
	for (i = 0; i < lattice_size; i++)
	{
		lattice[i] = 0.35;
	}

	path_integral(lattice, steps, 0.09, ho, cor);

    cor.analyze_data_sets();

	std::cout << cor;
}
