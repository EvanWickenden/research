
#include "lattice.h"

std::ostream& operator<<(std::ostream& s, const Lattice& lattice)
{
	int i;
	s << "[";
	for (i = 0; i < lattice.N; i++)
	{
		s << lattice[i] << ", ";
	}
	s << "]" << std::endl;
	return s;
}
