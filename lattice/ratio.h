#ifndef __RATIO_H__
#define __RATIO_H__


struct Ratio
{
	long numerator;
	long denominator;

	Ratio(long numerator, long denominator);

	Ratio& operator++(); 	/* prefix; increment numeratorerator */
	Ratio& operator++(int); /* postfix; increment denominatorominator */

	float operator()();
};


#endif
