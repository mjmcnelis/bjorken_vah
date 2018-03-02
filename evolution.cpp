
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "evolution.hpp"
#include "qcd.hpp"
#include "anisotropic_transport.hpp"




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                         EVOLUTION EQUATIONS                      ::
//                                                                  ::
//     Derivatives of hydrodynamic quantities with respect to       ::
//	   longitudinal proper time for bjorken flow. 		            ::
//                                                                  ::
//        dTtt_dtau	   dTtx_dtau	dTty_dtau    dTtn_dtau          ::
//																	::
//		  dpl_dtau     dpt_dtau                                     ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double dTtt_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	//double T = effectiveTemperature(e);
	//double t2 = t*t;
	//double z3 = ut/t/sqrt(1.0 + ux*ux + uy*uy);

	//double pt_macro = p + pt - equilibriumKineticPressure(T);
	//double pl_macro = p + pl - equilibriumKineticPressure(T);
	//double mbar = z_Quasiparticle(T);
	//double pkinetic = I21_function(T,mbar);
	//double Beq = equilibriumBquasi(T);

	//double pt_macro = pt - B;
	//double pl_macro = pl - B;

	//double Lnn = (pl-pt)*z3*z3;

	//double Tnn = (e+pt_macro)*un*un + pt_macro/t2 + Lnn;

	//return - (Ttt + t2*Tnn) / t;    // or - (Ttt + pl_macro) / t

	return -(Ttt + pl)/t;
}


double dTtx_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - Ttx / t;
}


double dTty_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - Tty / t;
}


double dTtn_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - 3.0 * Ttn / t;
}


double dpl_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	// thermodynamic functions
	double T = effectiveTemperature(e); 		         // temperature
	double s = (e+p)/T;									 // entropy density
	double etas = shearViscosityToEntropyDensity(T);     // eta / s
	double zetas = bulkViscosityToEntropyDensity(T);     // zeta / s
	double betapi = beta_shear(T);						 // eta / tau_pi
	double betabulk = beta_bulk(T);						 // zeta / tau_Pi
	double m = T * z_Quasiparticle(T);				     // m(T)
	double mbar = m/lambda;       					     // m(T) / lambda

	// relaxation rates
	double taupiInv = betapi/s/etas;
	double taubulkInv = betabulk/s/zetas;

	// kinetic pressure
	double pavg = (2.0*pt + pl) / 3.0;			         // average kinetic pressure
	double dp = pl - pt;						         // pressure anisotropy
	//double peq = equilibriumKineticPressure(T);

	// anisotropic functions
	double I240 = I240_function(lambda,ax,az,mbar);
	double I020 = I020_function(lambda,ax,az,mbar);

	// transport coefficient
	//double zetaLL = I240 - 3.0*pl + I020*(e+pl-B)*mdmde_Quasiparticle(e);
	double zetaLL = I240 - 3.0*(pl+b) + mdmde_Quasiparticle(e)*(e+pl)*(I020 + (2.0*pt+pl-e+4.0*b)/(m*m));

	// pldot
	//return taubulkInv*(peq-pavg) - 2.0*taupiInv*dp/3.0 + zetaLL/t;
	return taubulkInv*(p-pavg) - 2.0*taupiInv*dp/3.0 + zetaLL/t;
}


// time derivative of kinetic transverse pressure
double dpt_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	double T = effectiveTemperature(e);
	double s = (e+p)/T;
	double etas = shearViscosityToEntropyDensity(T);
	double zetas = bulkViscosityToEntropyDensity(T);
	double betapi = beta_shear(T);
	double betabulk = beta_bulk(T);

	double m = T * z_Quasiparticle(T);
	double mbar = m/lambda;

	double taupiInv = betapi/s/etas;
	double taubulkInv = betabulk/s/zetas;

	double pavg = (2.0*pt + pl) / 3.0;
	double dp = pl - pt;
	//double peq = equilibriumKineticPressure(T);

	double I221 = I221_function(lambda,ax,az,mbar);
	double I001 = I001_function(lambda,ax,az,mbar);

	//double zetaTL = I221 - pt + I001*(e+pl-B)*mdmde_Quasiparticle(e);
	double zetaTL = I221 - (pt+b) + mdmde_Quasiparticle(e)*(e+pl)*(I001 + (2.0*pt+pl-e+4.0*b)/(m*m));

	//return taubulkInv*(peq-pavg) + taupiInv*dp/3.0 + zetaTL/t;
	return taubulkInv*(p-pavg) + taupiInv*dp/3.0 + zetaTL/t;
}


// time derivative of non-equilibrium component of mean field dB
double db_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	double T = effectiveTemperature(e);
	double cs2 = speedOfSoundSquared(e);
	double s = (e+p)/T;
	double zetas = bulkViscosityToEntropyDensity(T);
	double betabulk = beta_bulk(T);

	double taubulkInv = betabulk/s/zetas;
	double m = T * z_Quasiparticle(T);

	double beq = equilibriumBquasi(T);

	// return -taubulkInv*(B-Beq) - (2.0*pt+pl-e+B)*mdmde_Quasiparticle(e)*(e+pl-B)/(t*m*m);
	return -taubulkInv*(b-beq) - (2.0*pt + pl - e + 4.0*b)*mdmde_Quasiparticle(e)*(e+pl)/(t*m*m);
}












