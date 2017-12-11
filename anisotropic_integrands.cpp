
#include <stdlib.h>
#include <math.h>
#include "anisotropic_integrands.hpp"
#include "rfunctions.hpp"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                  ANISOTROPIC FUNCTION INTEGRANDS                 ::
//                                                                  ::
//     Momentum bar integrand of the anisotropic functions.         ::
//	   (written in gla form) Integrate with Gauss_Aniso_1D.         ::
//     (prefactors pulled out)			                			::
//																	::
//     I240_integrand                           I221_integrand      ::
//                                                                  ::
//     Ea_integrand        PTa_integrand        PLa_integrand       ::
//     I2001_integrand     I2011_integrand      I2201_integrand     ::
//     I401m1_integrand    I420m1_integrand                         ::
//     I402m1_integrand    I421m1_integrand     I440m1_integrand    ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double I240_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R240(pbar,ax,az,mbar) * exp(pbar-sqrt(pbar*pbar + mbar*mbar));
}


double I221_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R221(pbar,ax,az,mbar) * exp(pbar-sqrt(pbar*pbar + mbar*mbar));
}


double I020_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R020(pbar,ax,az,mbar) * exp(pbar-sqrt(pbar*pbar + mbar*mbar));
}


double I001_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R001(pbar,ax,az,mbar) * exp(pbar-sqrt(pbar*pbar + mbar*mbar));
}


double Ea_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R200(pbar,ax,az,mbar) * exp(pbar - sqrt(pbar*pbar + mbar*mbar));
}


double PTa_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R201(pbar,ax,az,mbar) * exp(pbar - sqrt(pbar*pbar + mbar*mbar));
}


double PLa_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 2)
	return pbar * R220(pbar,ax,az,mbar) * exp(pbar - sqrt(pbar*pbar + mbar*mbar));
}


double I2001_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	return Eabar * R200(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


double I2011_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	return Eabar * R201(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


double I2201_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double Eabar = sqrt(pbar*pbar + mbar*mbar);
	return Eabar * R220(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


// double I401m1_integrand(double pbar, double ax, double az, double mbar)
// {
// 	// gauss laguerre (a = 3)
// 	double pbar2 = pbar * pbar;
// 	double Eabar = sqrt(pbar2 + mbar*mbar);
// 	return (pbar2 / Eabar) * R401(pbar,ax,az,mbar) * exp(pbar - Eabar);
// }


// double I420m1_integrand(double pbar, double ax, double az, double mbar)
// {
// 	// gauss laguerre (a = 3)
// 	double pbar2 = pbar * pbar;
// 	double Eabar = sqrt(pbar2 + mbar*mbar);
// 	return (pbar2 / Eabar) * R420(pbar,ax,az,mbar) * exp(pbar - Eabar);
// }


double I402m1_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double pbar2 = pbar * pbar;
	double Eabar = sqrt(pbar2 + mbar*mbar);
	return (pbar2 / Eabar) * R402(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


double I421m1_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double pbar2 = pbar * pbar;
	double Eabar = sqrt(pbar2 + mbar*mbar);
	return (pbar2 / Eabar) * R421(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


double I440m1_integrand(double pbar, double ax, double az, double mbar)
{
	// gauss laguerre (a = 3)
	double pbar2 = pbar * pbar;
	double Eabar = sqrt(pbar2 + mbar*mbar);
	return (pbar2 / Eabar) * R440(pbar,ax,az,mbar) * exp(pbar - Eabar);
}


double I21_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 2) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2 + mbar*mbar);

	return (pbar*pbar)/(Ebar) * exp(pbar-Ebar);
}


double I20_integrand(double pbar, double mbar)
{
	// gauss laguerre (a = 2) boltzmann
	double pbar2 = pbar * pbar;
	double Ebar = sqrt(pbar2 + mbar*mbar);

	return Ebar * exp(pbar-Ebar);
}










































