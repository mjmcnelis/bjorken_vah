
#include <stdlib.h>
#include "gauss_integration.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                         1D GAUSS INTEGRATION                       ::
//                                                                    ::
//     Compute 1D anisotropic integrals over radial momentum          ::
//     bar using Gauss Laguerre quadrature.                           ::
//                                                                    ::
//                           Gauss_Aniso_1D                           ::
//                                                                    ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double Gauss_Aniso_1D(double aniso_1D_integrand(double pbar, double ax, double az, double mbar), double * pbar_root, double * pbar_weight, const int pbar_pts, double ax, double az, double mbar)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * aniso_1D_integrand(pbar_root[k], ax, az, mbar);
	return sum;
}

double Gauss_Thermal_1D(double thermal_1D_integrand(double pbar, double mbar), double * pbar_root, double * pbar_weight, int pbar_pts, double mbar)
{
	double sum = 0.0;
	for(int k = 0; k < pbar_pts; k++) sum += pbar_weight[k] * thermal_1D_integrand(pbar_root[k], mbar);
	return sum;
}