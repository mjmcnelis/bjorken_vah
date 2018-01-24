
#include <stdlib.h>

#ifndef ANISOTROPICVARIABLES_H

#define ANISOTROPICVARIABLES_H


typedef enum {newton, broyden} jacobian;

void get_anisotropic_variables(double e, double pl, double pt, double B, double *lambda, double *ax, double *az);


#endif