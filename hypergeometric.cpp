
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hypergeometric.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                      HYPERGEOMETRIC FUNCTIONS                    ::
//                                                                  ::
//      Hypergeometric functions that appear in the R-functions     ::
//      of the anisotropic integrands (Ea,I2001), (PTa,I2011),      ::
//      (PLa,I2201), I401-1, I420-1, I402-1, I421-1, I441-1,        ::
//		respectively                                                ::
//                                                                  ::
//      tfuncE     tfuncPT     tfuncPL     tfunc401     tfunc420    ::
//                                                                  ::
//      tfunc402   tfunc421    tfunc440    tfunc240     tfunc221    ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double tfuncE(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = 1.0 + (1.0 + x) * atan(sqrt(x))/sqrt(x);
	else if(x < 0.0 && x > -1.0)
		result = 1.0 + (1.0 + x) * atanh(sqrt(-x))/sqrt(-x);
	else if(x == 0.0)
		result = 2.0;
	else if(x <= -1.0)
		throw "tfuncE outside domain";
	return result;
}

double tfuncPT(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// can I get rounding error propagation when |x| is very small?
		// yes, accuracy to 4/3 turns around when |x| = 1e-9 and lower
		result = (1.0 + (x - 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0. && x > -1.0)
		result = (1.0 + (x - 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 4.0 / 3.0;
	else if(x <= -1.0)
		throw "tfuncPT outside domain";
	return result;
}

double tfuncPL(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// rounding error propagation when |x| is very small
		// accuracy to 2/3 turns around when |x| = 1e-9 and lower
		result = (-1.0 + (x + 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0. && x > -1.)
		result = (-1.0 + (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 2.0 / 3.0;
	else if(x <= -1.0)
		throw "tfuncPL outside domain";
	return result;
}

double tfunc401(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// accuracy to 4/3 turns around when |x| = 1e-10 and lower (rounding error)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < 0.0 && x > -1.0)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x == 0.0)
		result = 4.0 / 3.0;
	else if(x <= -1.0)
		throw "tfunc401 outside domain";
	return result;
}

double tfunc420(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < 0.0 && x > -1.0)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x == 0.0)
		result = 2.0 / 3.0;
	else if(x <= -1.0)
		throw "tfunc420 outside domain";
	return result;
}

double tfunc402(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 16/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 16.0 / 15.0;
	else if(x <= -1.0)
		throw "tfunc402 outside domain";
	return result;
}

double tfunc421(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 4/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 4.0 / 15.0;
	else if(x <= -1.0)
		throw "tfunc421 outside domain";
	return result;
}

double tfunc440(double x)
{
	double result = 0.0;
	if(x > 0.0)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 2/5 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < 0.0 && x > -1.0)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x == 0.0)
		result = 0.4;
	else if(x <= -1.0)
		throw "tfunc440 outside domain";
	return result;
}


double tfunc240(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = (3.0+2.0*x - 3.0*(1.0+x)*atan(sqrt(x))/sqrt(x)) / (x*x);
	else if(x < 0.0 && x > -1.0)
		result = (3.0+2.0*x - 3.0*(1.0+x)*atanh(sqrt(-x))/sqrt(-x)) / (x*x);
	else if(x == 0.0)
		result = 0.4;
	else if(x <= -1.0)
		throw "tfunc240 outside domain";
	return result;
}


double tfunc221(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = (-3.0 + (3.0+x)*atan(sqrt(x))/sqrt(x)) / (x*x);
	else if(x < 0.0 && x > -1.0)
		result = (-3.0 + (3.0+x)*atanh(sqrt(-x))/sqrt(-x)) / (x*x); // fixed bug on 12/8
	else if(x == 0.0)
		result = 4.0/15.0;
	else if(x <= -1.0)
		throw "tfunc221 outside domain";
	return result;
}

double tfunc020(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = 2.0 * (1.0 - atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0.0 && x > -1.0)
		result = 2.0 * (1.0 - atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 2.0/3.0;
	else if(x <= -1.0)
		throw "tfunc020 outside domain";
	return result;
}

double tfunc001(double x)
{
	double result = 0.0;
	if(x > 0.0)
		result = 2.0 * (-1.0/(1.0+x) + atan(sqrt(x))/sqrt(x)) / x;
	else if(x < 0.0 && x > -1.0)
		result = 2.0 * (-1.0/(1.0+x) + atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x == 0.0)
		result = 4.0/3.0;
	else if(x <= -1.0)
		throw "tfunc001 outside domain";
	return result;
}
