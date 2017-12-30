
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;
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
	double a = 0.01;

	if(x > a)
		result = 1.0 + (1.0 + x) * atan(sqrt(x))/sqrt(x);
	else if(x < -a && x > -1.0)
		result = 1.0 + (1.0 + x) * atanh(sqrt(-x))/sqrt(-x);
	else if(x >= -a && x <= a)
		result = 2.0 + x*(0.6666666666666667 + x*(-0.1333333333333333 +
			x*(0.05714285714285716 + x*(-0.031746031746031744 + x*(0.020202020202020193 +
			x*(-0.013986013986013984 + (0.010256410256410262 - 0.00784313725490196*x)*x))))));
	else if(x <= -1.0)
		throw "tfuncE(x) is out of bounds";
	return result;
}

double tfuncPT(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// can I get rounding error propagation when |x| is very small?
		// yes, accuracy to 4/3 turns around when |x| = 1e-9 and lower
		result = (1.0 + (x - 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < -a && x > -1.0)
		result = (1.0 + (x - 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x >= -a && x <= a)
		result = 1.3333333333333333 + x*(-0.5333333333333333 + x*(0.34285714285714286 +
         	x*(-0.25396825396825395 + x*(0.20202020202020202 + x*(-0.16783216783216784 +
         	x*(0.14358974358974358 + (-0.12549019607843137 + 0.11145510835913312*x)*x))))));
	else if(x <= -1.0)
		throw "tfuncPT(x) is out of bounds";

	return result;
}

double tfuncPL(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// rounding error propagation when |x| is very small
		// accuracy to 2/3 turns around when |x| = 1e-9 and lower
		result = (-1.0 + (x + 1.0) * atan(sqrt(x))/sqrt(x)) / x;
	else if(x < -a && x > -1.0)
		result = (-1.0 + (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x >= -a && x <= a)
		result = 0.6666666666666667 + x*(-0.1333333333333333 + x*(0.05714285714285716 +
			x*(-0.031746031746031744 + x*(0.020202020202020193 + x*(-0.013986013986013984 +
			x*(0.010256410256410262 + (-0.00784313725490196 + 0.006191950464396287*x)*x))))));
	else if(x <= -1.0)
		throw "tfuncPL(x) is out of bounds";

	return result;
}

double tfunc401(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// accuracy to 4/3 turns around when |x| = 1e-10 and lower (rounding error)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < -a && x > -1.0)
		result = (1.0 + 3.0*x + (3.0*x - 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x >= -a && x <= a)
		result = 1.3333333333333335 + x*(0.5333333333333333 + x*(-0.11428571428571427 +
         	x*(0.05079365079365082 + x*(-0.02886002886002885 + x*(0.01864801864801864 +
         	x*(-0.01305361305361305 + (0.00965309200603319 - 0.007430340557275546*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc401(x) is out of bounds";

	return result;
}

double tfunc420(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x);
	else if(x < -a && x > -1.0)
		result = (x - 1.0 + (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x);
	else if(x >= -a && x <= a)
		result = 0.6666666666666667 + x*(0.13333333333333336 + x*(-0.019047619047619035 +
         x*(0.006349206349206354 + x*(-0.0028860028860028877 + x*(0.0015540015540015523 +
         x*(-0.0009324009324009307 + (0.0006033182503770752 - 0.00041279669762641844*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc420(x) is out of bounds";

	return result;
}

double tfunc402(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 16/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < -a && x > -1.0)
		result = (3.0 * (x - 1.0) + (3.0 + (x * (3.0*x - 2.0))) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x >= -a && x <= a)
		result = 1.0666666666666667 + x*(-0.45714285714285713 + x*(0.3047619047619048 +
         x*(-0.23088023088023085 + x*(0.1864801864801865 + x*(-0.15664335664335666 +
         x*(0.13514328808446457 + (-0.11888544891640866 + 0.10614772224679345*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc402(x) is out of bounds";

	return result;
}

double tfunc421(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 4/15 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < -a && x > -1.0)
		result = (3.0 + x + (x + 1.0) * (x - 3.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x >= -a && x <= a)
		result = 0.2666666666666666 + x*(-0.0761904761904762 + x*(0.0380952380952381 +
			x*(-0.023088023088023088 + x*(0.015540015540015537 + x*(-0.011188811188811189 +
            x*(0.00844645550527904 + (-0.006604747162022705 + 0.005307386112339673*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc421(x) is out of bounds";

	return result;
}

double tfunc440(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		// more prone to rounding error when x is small due to x^2 / x^2
		// accuracy to 2/5 turns around when |x| = 1e-6 and lower  (can be numerically unstable)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atan(sqrt(x))/sqrt(x)) / (4.0*x*x);
	else if(x < -a && x > -1.0)
		result = (-(3.0 + 5.0*x) + 3.0 * (x + 1.0) * (x + 1.0) * atanh(sqrt(-x))/sqrt(-x)) / (4.0*x*x);
	else if(x >= -a && x <= a)
		result = 0.4 + x*(-0.057142857142857106 + x*(0.019047619047619063 +
			x*(-0.008658008658008663 + x*(0.004662004662004657 + x*(-0.002797202797202792 +
            x*(0.0018099547511312257 + (-0.0012383900928792553 + 0.0008845643520566139*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc440(x) is out of bounds";

	return result;
}


double tfunc240(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = (3.0+2.0*x - 3.0*(1.0+x)*atan(sqrt(x))/sqrt(x)) / (x*x);
	else if(x < -a && x > -1.0)
		result = (3.0+2.0*x - 3.0*(1.0+x)*atanh(sqrt(-x))/sqrt(-x)) / (x*x);
	else if(x >= -a && x <= a)
		result = 0.4 + x*(-0.17142857142857149 + x*(0.09523809523809523 +
			x*(-0.06060606060606058 + x*(0.04195804195804195 + x*(-0.030769230769230785 +
            x*(0.023529411764705882 + (-0.01857585139318886 + 0.015037593984962405*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc240(x) is out of bounds";

	return result;
}


double tfunc221(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = (-3.0 + (3.0+x)*atan(sqrt(x))/sqrt(x)) / (x*x);
	else if(x < -a && x > -1.0)
		result = (-3.0 + (3.0+x)*atanh(sqrt(-x))/sqrt(-x)) / (x*x); // fixed bug on 12/8
	else if(x >= -a && x <= a)
		result = 0.2666666666666668 + x*(-0.22857142857142854 + x*(0.19047619047619047 +
			x*(-0.1616161616161616 +  x*(0.13986013986013987 + x*(-0.12307692307692308 +
            x*(0.10980392156862744 + (-0.09907120743034055 + 0.09022556390977443*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc221(x) is out of bounds";

	return result;
}

double tfunc020(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = 2.0 * (1.0 - atan(sqrt(x))/sqrt(x)) / x;
	else if(x < -a && x > -1.0)
		result = 2.0 * (1.0 - atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x >= -a && x <= a)
		result = 0.6666666666666666 + x*(-0.4 + x*(0.2857142857142857 +
         x*(-0.2222222222222222 + x*(0.18181818181818182 + x*(-0.15384615384615385 +
         x*(0.13333333333333333 + (-0.11764705882352941 + 0.10526315789473684*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc020(x) is out of bounds";

	return result;
}

double tfunc001(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = 2.0 * (-1.0/(1.0+x) + atan(sqrt(x))/sqrt(x)) / x;
	else if(x < -a && x > -1.0)
		result = 2.0 * (-1.0/(1.0+x) + atanh(sqrt(-x))/sqrt(-x)) / x;
	else if(x >= -a && x <= a)
		result = 1.3333333333333335 + x*(-1.6 + x*(1.7142857142857144 +
         x*(-1.7777777777777777 + x*(1.8181818181818181 + x*(-1.8461538461538463 +
         x*(1.8666666666666667 + x*(-1.8823529411764706 + 1.8947368421052633*x)))))));
	else if(x <= -1.0)
		throw "tfunc001(x) is out of bounds";

	return result;
}

double tfunc000(double x)
{
	double result = 0.0;
	double a = 0.01;

	if(x > a)
		result = 2.0 * atan(sqrt(x))/sqrt(x);
	else if(x < -a && x > -1.0)
		result = 2.0 * atanh(sqrt(-x))/sqrt(-x);
	else if(x >= -a && x <= a)
		result = 2.0 + x*(-0.6666666666666666 + x*(0.4 + x*(-0.2857142857142857 +
            x*(0.2222222222222222 + x*(-0.18181818181818182 +
                  x*(0.15384615384615385 + (-0.13333333333333333 + 0.11764705882352941*x)*x))))));
	else if(x <= -1.0)
		throw "tfunc000(x) is out of bounds";

	return result;
}
