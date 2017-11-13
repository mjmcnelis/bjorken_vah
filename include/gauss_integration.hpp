
#include <stdlib.h>

#ifndef GAUSS_INTEGRATION_H

#define GAUSS_INTEGRATION_H


double Gauss_Aniso_1D(double aniso_1D_integrand(double pbar, double ax, double az, double mbar), double * pbar_root, double * pbar_weight, const int pbar_pts, double ax, double az, double mbar);


#endif
