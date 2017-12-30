
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
	double t2 = t*t;
	double z3 = ut/t/sqrt(1.0 + ux*ux + uy*uy);

	//double pt_macro = p + pt - equilibriumKineticPressure(T);
	//double pl_macro = p + pl - equilibriumKineticPressure(T);
	//double mbar = z_Quasiparticle(T);
	//double pkinetic = I21_function(T,mbar);
	//double Beq = equilibriumBquasi(T);
	double pt_macro = pt - B;
	double pl_macro = pl - B;

	double Lnn = (pl-pt)*z3*z3;

	double Tnn = (e+pt_macro)*un*un + pt_macro/t2 + Lnn;

	double Tttdot = - (Ttt + t2*Tnn) / t;    // or - (Ttt + pl_macro) / tau

	return Tttdot;
}


double dTtx_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - Ttx / t;
}


double dTty_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - Tty / t;
}


double dTtn_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	return - 3.0 * Ttn / t;
}


// time derivative of kinetic longitudinal pressure
double dpl_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	// lattice qcd
	double T = effectiveTemperature(e); 		         // temperature
	double s = (e+p)/T;									 // entropy density
	double etas = shearViscosityToEntropyDensity(T);     // eta / s
	double zetas = bulkViscosityToEntropyDensity(T);     // zeta / s
	double betapi = beta_shear(T);						 // eta / tau_pi
	double betaPi = beta_bulk(T);						 // zeta / tau_Pi
	double mbar = z_Quasiparticle(T) * (T/lambda);       // m(T) / lambda


	// relaxation rates
	double taupiInv = betapi/s/etas;                     // shear relaxation rate
	double tauPiInv = betaPi/s/zetas;                    // bulk relaxation rate

	//double cs2 = speedOfSoundSquared(e);
	//double b2 = (1.0/3.0 - cs2);
	//double taupiInv = 0.2 * T / etas;
	//double tauPiInv = 15.0 * b2 * b2 * T / zetas;


	double pavg = (2.0*pt + pl) / 3.0;			         // average kinetic pressure
	double dp = pl - pt;						         // pressure anisotropy
	double mbar_eq = z_Quasiparticle(T);                 // m(T) / T
	double peq = I21_function(T,mbar_eq);                // equilibrium kinetic pressure

	double I240 = I240_function(lambda,ax,az,mbar);
	double I020 = I020_function(lambda,ax,az,mbar);

	double edot = dTtt_dt(Ttt, Ttx, Tty, Ttn, pl, pt, B, ut, ux, uy, un, e, p, lambda, ax, az, t);


	//double dTde = 1.0 / derivativeEnergyDensityWithRespectToTemperature(T);
	//double quasi_term = I020 * mdmdT_Quasiparticle(T) * dTde * edot;

	double quasi_term = I020 * mdmde_Quasiparticle(e) * edot;

	// set edot = (e_n - e_n-1) / dtau backward difference; feed in e, ep and dtau;

	// kinetic pl relaxation equation
	double pldot = tauPiInv*(peq-pavg) - 2.0*taupiInv*dp/3.0 + (I240-3.0*pl)/t - quasi_term;

	return pldot;
}


// time derivative of kinetic transverse pressure
double dpt_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	double T = effectiveTemperature(e);
	double s = (e+p)/T;
	double etas = shearViscosityToEntropyDensity(T);
	double zetas = bulkViscosityToEntropyDensity(T);
	double betapi = beta_shear(T);
	double betaPi = beta_bulk(T);
	double mbar = z_Quasiparticle(T) * (T/lambda);

	double taupiInv = betapi/s/etas;
	double tauPiInv = betaPi/s/zetas;

	//double cs2 = speedOfSoundSquared(e);
	//double b2 = (1.0/3.0 - cs2);
	//double taupiInv = 0.2 * T / etas;
	//double tauPiInv = 15.0 * b2 * b2 * T / zetas;

	double pavg = (2.0*pt + pl) / 3.0;
	double dp = pl - pt;
	double mbar_eq = z_Quasiparticle(T);
	double peq = I21_function(T,mbar_eq);

	double I221 = I221_function(lambda,ax,az,mbar);
	double I001 = I001_function(lambda,ax,az,mbar);

	double edot = dTtt_dt(Ttt, Ttx, Tty, Ttn, pl, pt, B, ut, ux, uy, un, e, p, lambda, ax, az, t);

	//double dTde = 1.0 / derivativeEnergyDensityWithRespectToTemperature(T);
	//double quasi_term = I001 * mdmdT_Quasiparticle(T) * dTde * edot;

	double quasi_term = I001 * mdmde_Quasiparticle(e) * edot;

	// kinetic pl relaxation equation
	double ptdot = tauPiInv*(peq-pavg) + taupiInv*dp/3.0 + (I221-pt)/t - quasi_term;

	return ptdot;
}


// time derivative of non-equilibrium component of mean field dB
double dB_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double B, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t)
{
	double T = effectiveTemperature(e);
	double s = (e+p)/T;
	double etas = shearViscosityToEntropyDensity(T);
	double zetas = bulkViscosityToEntropyDensity(T);
	double betapi = beta_shear(T);
	double betaPi = beta_bulk(T);
	double mbar = z_Quasiparticle(T) * (T/lambda);

	double taupiInv = betapi/s/etas;
	double tauPiInv = betaPi/s/zetas;
	double mbar_eq = z_Quasiparticle(T);
	double peq = I21_function(T,mbar_eq);

	double I000 = I000_function(lambda,ax,az,mbar);
	double I00_eq = I00_function(T,mbar_eq);

	double edot = dTtt_dt(Ttt, Ttx, Tty, Ttn, pl, pt, B, ut, ux, uy, un, e, p, lambda, ax, az, t);


	double Beq = equilibriumBquasi(T);
	//double dB = B - Beq;
	//double Eeq = I20_function(T,mbar_eq);
	double m_eq = mbar_eq * T;
	//double quasi_term = -(I000 - I00_eq) * mdmde_Quasiparticle(e) * edot;
	//double quasi_term = - I000 * mdmde_Quasiparticle(e) * edot;
	//double quasi_term = (2.0*pt + pl - Eeq + dB) * mdmde_Quasiparticle(e) * edot / (m_eq * m_eq);
	double quasi_term = (2.0*pt + pl - e + B) * mdmde_Quasiparticle(e) * edot / (m_eq * m_eq);

	//double quasi_term = (2.0*pt + pl - 3.0*peq + dB) * mdmdT_Quasiparticle(T) * dTde * edot / (mbar_eq * mbar_eq * T * T);

	// kinetic pl relaxation equation
	double Bdot = -tauPiInv*(B-Beq) + quasi_term;

	return Bdot;
}














