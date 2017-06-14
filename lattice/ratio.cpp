
#include "ratio.h"

#include "log.h"

Ratio::Ratio(long numerator, long denominator) :
	numerator(numerator),
	denominator(denominator)
{}

/* prefix */
Ratio& Ratio::operator++()
{
	numerator++;
	return *this;
}

Ratio& Ratio::operator++(int)
{
	denominator++;
	return *this;
}

float Ratio::operator()()
{
	return numerator / (float) denominator;
}
