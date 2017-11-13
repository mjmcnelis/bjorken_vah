
#include <stdlib.h>

#ifndef ANISOTROPIC_INTEGRANDS_H

#define ANISOTROPIC_INTEGRANDS_H


double I240_integrand(double pbar, double ax, double az, double mbar);

double I221_integrand(double pbar, double ax, double az, double mbar);


double Ea_integrand(double pbar, double ax, double az, double mbar);
double I2001_integrand(double pbar, double ax, double az, double mbar);

double PTa_integrand(double pbar, double ax, double az, double mbar);
double I2011_integrand(double pbar, double ax, double az, double mbar);

double PLa_integrand(double pbar, double ax, double az, double mbar);
double I2201_integrand(double pbar, double ax, double az, double mbar);

double I401m1_integrand(double pbar, double ax, double az, double mbar);
double I420m1_integrand(double pbar, double ax, double az, double mbar);

double I402m1_integrand(double pbar, double ax, double az, double mbar);
double I421m1_integrand(double pbar, double ax, double az, double mbar);
double I440m1_integrand(double pbar, double ax, double az, double mbar);


#endif