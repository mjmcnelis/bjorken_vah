
#include <stdlib.h>
#include <math.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "evolution.hpp"
#include "qcd.hpp"
#include "anisotropic_functions.hpp"
#include "anisotropic_integrands.hpp"
#include "gauss_integration.hpp"



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                   ANISOTROPIC TRANSPORT FUNCTIONS                ::
//                                                                  ::
//     Moments of anisotropic distribution function that appear     ::
//     in transport coefficients of anisotropic hydrodynamics.      ::
//     													            ::
//       				  I240 			I221                        ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double I240_function(double lambda, double ax, double az, double mbar)
{
		const int pbar_pts = 16;
		// a = 2 gauss laguerre roots/weight
		double pbar_root[pbar_pts] = {0.377613508344741,1.01749195760257,1.94775802042424,3.17692724488987,4.7162400697918,6.58058826577491,8.78946527064707,11.3683230828333,14.350626727437,17.7810957248416,21.7210847965713,26.2581386751111,31.5245960042758,37.7389210025289,45.3185461100898,55.3325835388358};
		double pbar_weight[pbar_pts] = {0.0486064094670787,0.29334739019044,0.583219363383551,0.581874148596173,0.33818053747379,0.12210596394498,0.0281146258006637,0.00414314919248226,0.000385648533767438,2.20158005631091e-05,7.34236243815652e-07,1.32646044204804e-08,1.15266648290843e-10,3.94706915124609e-13,3.63797825636053e-16,3.45457612313612e-20};

		double g = 51.4103536012791;

		double factor_I240 = g * (ax*ax) * (az*az*az*az*az) * (lambda*lambda*lambda*lambda) / (4.0*M_PI*M_PI);

		double answer = factor_I240 * Gauss_Aniso_1D(I240_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		// check conformal formula later:

		return answer;
}



double I221_function(double lambda, double ax, double az, double mbar)
{
		const int pbar_pts = 16;
		// a = 2 gauss laguerre roots/weights
		double pbar_root[pbar_pts] = {0.377613508344741,1.01749195760257,1.94775802042424,3.17692724488987,4.7162400697918,6.58058826577491,8.78946527064707,11.3683230828333,14.350626727437,17.7810957248416,21.7210847965713,26.2581386751111,31.5245960042758,37.7389210025289,45.3185461100898,55.3325835388358};
		double pbar_weight[pbar_pts] = {0.0486064094670787,0.29334739019044,0.583219363383551,0.581874148596173,0.33818053747379,0.12210596394498,0.0281146258006637,0.00414314919248226,0.000385648533767438,2.20158005631091e-05,7.34236243815652e-07,1.32646044204804e-08,1.15266648290843e-10,3.94706915124609e-13,3.63797825636053e-16,3.45457612313612e-20};

		const double g = 51.4103536012791; // degeneracy factor g (nf = 3 flavors)

		double factor_I221 = g * (ax*ax*ax*ax) * (az*az*az) * (lambda*lambda*lambda*lambda) / (8.0*M_PI*M_PI);;

		double answer = factor_I221 * Gauss_Aniso_1D(I221_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		// check conformal formula later:

		return answer;
}









