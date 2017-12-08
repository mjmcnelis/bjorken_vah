
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "evolution.hpp"
#include "qcd.hpp"
#include "anisotropic_functions.hpp"




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


double dTtt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	double T = effectiveTemperature(e);
	double tau2 = tau*tau;
	double z3 = ut/tau/sqrt(1.0 + ux*ux + uy*uy);

	double pt_macro = p + pt - equilibriumKineticPressure(T);
	double pl_macro = p + pl - equilibriumKineticPressure(T);

	double Lnn = (pl-pt)*z3*z3;

	double Tnn = (e+pt_macro)*un*un + pt_macro/tau2 + Lnn;

	double Tttdot = - (Ttt + tau2*Tnn) / tau;    // or - (Ttt + pl_macro) / tau

	return Tttdot;
}


double dTtx_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	return - Ttx / tau;
}


double dTty_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	return - Tty / tau;
}


double dTtn_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	return - 3.0 * Ttn / tau;
}


// time derivative of kinetic longitudinal pressure
double dpl_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	// lattice qcd
	double T = effectiveTemperature(e); 		         // temperature
	double s = (e+p)/T;									 // entropy density
	double etas = shearViscosityToEntropyDensity(T);     // specific shear viscosity
	double zetas = bulkViscosityToEntropyDensity(T);     // specific bulk viscosity
	double betapi = beta_pi(T);						     // eta / tau_pi
	double betaPi = beta_Pi(T);						     // zeta / tau_Pi
	double mbar = z_Quasiparticle(T) * (T/lambda);       // m(T) / lambda


	// relaxation rates
	double taupiInv = betapi/s/etas;                     // shear relaxation rate
	double tauPiInv = betaPi/s/zetas;                    // bulk relaxation rate


	double pavg = (2.0*pt + pl) / 3.0;			         // average kinetic pressure
	double dp = pl - pt;						         // pressure anisotropy
	double peq = equilibriumKineticPressure(T);          // kinetic contribution only (no B(T))

	double I240 = I240_function(lambda,ax,az,mbar);


	// take out quasi term for now; if put back add to r.h.s: pldot += - quasi_term

	//double I020 = I020(lambda,ax,az,T);
	//double quasi_term = I020 * mdmdT_Quasiparticle(T); // quasiparticle contribution

	// shoot the quasi_term is not correct, I need T_dot or e_dot

	// requires m * dm/de * edot
	// rewrite mdmdT to mdmde
	// set edot = (e_n - e_n-1) / dtau backward difference; feed in e, ep and dtau;

	// kinetic pl relaxation equation
	double pldot = tauPiInv*(peq-pavg) - 2.0*taupiInv*dp/3.0 + (I240-3.0*pl)/tau;

	return pldot;
}


// time derivative of kinetic transverse pressure
double dpt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau)
{
	double T = effectiveTemperature(e);
	double s = (e+p)/T;
	double etas = shearViscosityToEntropyDensity(T);
	double zetas = bulkViscosityToEntropyDensity(T);
	double betapi = beta_pi(T);
	double betaPi = beta_Pi(T);
	double mbar = z_Quasiparticle(T) * (T/lambda);

	double taupiInv = betapi/s/etas;
	double tauPiInv = betaPi/s/zetas;

	double pavg = (2.0*pt + pl) / 3.0;
	double dp = pl - pt;
	double peq = equilibriumKineticPressure(T);

	double I221 = I221_function(lambda, ax, az, mbar);


	// take out quasi term for now; if put back add to r.h.s: pldot += - quasi_term

	//double I001 = I001(lambda,ax,az,T);
	//double quasi_term = I001 * mdmdT_Quasiparticle(T);   // quasiparticle contribution

	// shoot the quasi_term is not correct, I need T_dot or e_dot

	// requires m * dm/de * edot
	// rewrite mdmdT to mdmde
	// set edot = (e_n - e_n-1) / dtau backward difference; feed in e, ep and dtau;


	// kinetic pl relaxation equation
	double ptdot = tauPiInv*(peq-pavg) + taupiInv*dp/3.0 + (I221-pt)/tau;

	return ptdot;
}

















