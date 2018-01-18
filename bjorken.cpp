#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "qcd.hpp"
#include <ctime>
#include "evolution.hpp"
#include "inferredvariables.hpp"
#include "anisotropicvariables.hpp"
#include "anisotropic_transport.hpp"

//#include <gsl/gsl_sf.h>


#define GEV_TO_INVERSE_FM 5.067731
//temporary
const int alpha = 21;
const int gla_pts = 32;
double root_gla[alpha][gla_pts];
double weight_gla[alpha][gla_pts];

int load_gauss_laguerre_data()
{
  FILE *fp;
  stringstream laguerre_roots_weights;
  laguerre_roots_weights << "gla_roots_weights_" << gla_pts << "_points.txt";
  if((fp = fopen(laguerre_roots_weights.str().c_str(), "r")) == NULL)
  {
     return 1;
  }
  for(int i = 0; i < alpha; i++)
  {
   for(int j = 0; j < gla_pts; j++)
   {
      if(fscanf(fp, "%i %lf %lf", &i, &root_gla[i][j], &weight_gla[i][j])!= 3)
      	{
        	printf("error loading roots/weights of Gauss-Laguerre Quadradture at %d %d\n", i, j);
    		return 1;
    	}
   }
  }
  fclose(fp);
  return 0;
}




int main()
{
	clock_t begin;
    double duration;
    begin = clock();

 //    int num_error;
	// // Load gauss laguerre roots-weights
	// if((num_error = load_gauss_laguerre_data()) != 0)
	// {
	// 	fprintf(stderr, "Error loading gauss data (%d)!\n", num_error);
	// 	return 1;
	// }

	// for(int k = 0; k < gla_pts; k++) cout << setprecision(15) << root_gla[0][k] << ",";
	// printf("\n\n\n");
	// for(int k = 0; k < gla_pts; k++) cout << setprecision(15) << weight_gla[0][k] << ",";

	// Bjorken flow

	// input parameters
	double T0 = 0.6 * GEV_TO_INVERSE_FM;  // initial temperature in fm^-1
	double t0 = 0.25;				      // initial time in fm
	double tf = 50.0;					  // final time in fm


	// initialize temperature
	double T = T0;


	double mbar_eq = z_Quasiparticle(T);

	//double pkinetic = equilibriumKineticPressure(T0);

	//double pkinetic = I21_function(T,mbar_eq);
	double Beq = equilibriumBquasi(T);


	// initial flow profile
	double ut = 1.0;
	double ux = 0.0;
	double uy = 0.0;
	double un = 0.0;


	// initial energy density and pressure  (units = [fm^-4])
	double e0 = equilibriumEnergyDensity(T0);
	double e = e0;
	double p0 = equilibriumPressure(e0);
	double p = p0;

	double s0 = (e0+p0)/T0;
	double zetas0 = bulkViscosityToEntropyDensity(T0);    // initial specific bulk viscosity
	double etas0 = shearViscosityToEntropyDensity(T0);

	double taupi = (s0*etas0) / beta_shear(T0);           // quasiparticle relaxation times
	double taubulk = (s0*zetas0) / beta_bulk(T0);

	double piNS =  4.0 * (e0+p0) / (3.0*T0*t0) * shearViscosityToEntropyDensity(T0);
	double bulkNS = - (e0+p0) / (t0*T0) * bulkViscosityToEntropyDensity(T0);


	// initial T^{\tau\mu} components (units = [fm^-4])
	double Ttt = e0;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttn = 0.0;

	//double dB2nd = -3.0*taubulk*mdmdT_Quasiparticle(T)/pow(z_Quasiparticle(T),2)*speedOfSoundSquared(e)*(2.0*pt/3.0+pl/3.0-B-p)/(t0*T);
	//double B = Beq + dB2nd;

	//double plmacro = 1.0 * p;
	//double ptmacro = 1.0 * p;

	//double pl = plmacro + B;
	//double pt = ptmacro + B;

	// microscopic initialization
	// double lambda = 1.0 * T;
	// double ax = 1.0;
	// double az = 1.0;
	// double mbar = z_Quasiparticle(T) * (T/lambda);
	// double Ea = Ea_function(lambda, ax, az, mbar);
	// double pt = PTa_function(lambda, ax, az, mbar);
	// double pl = PLa_function(lambda, ax, az, mbar);
	// double B = e - Ea;  // perhaps this can be updated along with lambda, ax, az
	// double dB2nd = -3.0*taubulk*mdmdT_Quasiparticle(T)/pow(z_Quasiparticle(T),2)*speedOfSoundSquared(e)*(2.0*pt/3.0+pl/3.0-B-p)/(t0*T);


	// macroscopic initialization

	// limit PL/PT ~ 0.15 for zero bulk pressure
	//double PL = 0.209 * p;
	//double PT = 1.3955 * p;

	// so maybe I can adjust the nonequilibrium mean field until I can find a solution

	// regulation-scale total B until a solution can be inverted
	// or scale the dB only? the regulation should depend on the degree of pressure anisotropy

// Glasma initial conditions
	//double PL = 0.014925 * e / 3.0;
	//double PT = 1.4925 * e / 3.0;

	double PLPTratio = 0.01;

	double PT = (3.0/(2.0+PLPTratio)) * e / 3.0;
	double PL = (3.0 - 6.0/(2.0+PLPTratio)) * e / 3.0;

	PT = p;
	PL = p; 


	// double dB2nd = -3.0*taubulk*mdmdT_Quasiparticle(T)/pow(z_Quasiparticle(T),2)*speedOfSoundSquared(e)*(2.0*PT/3.0+PL/3.0-p)/(t0*T);
	// double B = 0.1141*(Beq + dB2nd);
	// cout << dB2nd << endl;

	double m = T * z_Quasiparticle(T);

	// double a = 3.0*taubulk*mdmde_Quasiparticle(e)*(e+PL)/(t0*m0*m0);
	// cout << 4.0/3.0*a << endl;

	double dBasy = -3.0*taubulk*mdmde_Quasiparticle(e)*(e+PL)*(2.0*PT/3.0+PL/3.0-p)/(t0*m*m) /
			(1.0 + 4.0*taubulk*mdmde_Quasiparticle(e)*(e+PL)/(t0*m*m));

	//double B = 0.142*(Beq + dBasy);



	double B = 1.0*(Beq + dBasy);

	//cout << (2.0*PT/3.0+PL/3.0-p) << endl;
	//cout << B - Beq << endl;
	//cout << dBasy << endl;
	double pl = PL + B;
	double pt = PT + B;
	double Ea = e - B;
	double lambda = 1.0 * T0;
	double ax = 1.0;
	double az = 1.0;



	try
	{
	get_anisotropic_variables(e, pl, pt, B, &lambda, &ax, &az);
	}
	catch (char const *excp)
	{
        	cout << "\nInitialization error: " << excp;
        	exit(-1);
    }

    double lambda0 = lambda;
    double ax0 = ax;
    double az0 = az;

    cout << "\n\n PL/PT = " << PL/PT << endl;

	cout << "\nlambda = " << lambda0 << ";" << endl;
	cout << "ax = " << ax0 << ";" << endl;
	cout << "az = " << az0 << ";" << endl;

	//exit(-1);
	//cout << p << endl;
	//cout << pkinetic - Beq << endl;


	// initialize mean field component (temporary)
	//double B = Beq;
	//double pl = I21_function(T0,mbar0);
	//double pt = I21_function(T0,mbar0);







	// intermediate and end values
	double Ttt_mid, Ttt_end;
	double Ttx_mid, Ttx_end;
	double Tty_mid, Tty_end;
	double Ttn_mid, Ttn_end;
	double pl_mid, pl_end;
	double pt_mid, pt_end;
	double B_mid, B_end;



	// time data (t = tau)
	double t = t0;
	const double dt = 0.005;
	const int n = floor((tf - t0) / dt);
	const int timesteps_per_write = 10;



	// Data files for plots
	ofstream eplot, piplot, bulkplot, plptplot;
	ofstream Tplot, lambdaplot, axplot, azplot;
	ofstream RpiInvplot, RbulkInvplot;
	ofstream piNSplot, bulkNSplot;
	ofstream taupiplot, taubulkplot;

	ofstream Bplot, Beqplot, dBasyplot;

	eplot.open("eplot_vah.dat", ios::out);
	piplot.open("piplot_vah.dat", ios::out);
	bulkplot.open("bulkplot_vah.dat", ios::out);
	plptplot.open("plptplot_vah.dat", ios::out);

	Tplot.open("Tplot_vah.dat", ios::out);
	lambdaplot.open("lambdaplot_vah.dat", ios::out);
	axplot.open("axplot_vah.dat", ios::out);
	azplot.open("azplot_vah.dat", ios::out);

	RpiInvplot.open("RpiInvplot_vah.dat", ios::out);
	RbulkInvplot.open("RbulkInvplot_vah.dat", ios::out);

	piNSplot.open("piNSplot_vah.dat", ios::out);
	bulkNSplot.open("bulkNSplot_vah.dat", ios::out);

	taupiplot.open("taupiplot_vah.dat", ios::out);
	taubulkplot.open("taubulkplot_vah.dat", ios::out);

	Bplot.open("Bplot_vah.dat", ios::out);
	Beqplot.open("Beqplot_vah.dat", ios::out);
	dBasyplot.open("dBasyplot_vah.dat", ios::out);


	eplot << "t [fm]" << "\t\t" << "e/e0" << endl << setprecision(5) << t << "\t\t" << e/e0 << endl;
	piplot << "t [fm]" << "\t\t" << "pi [fm^-4]" << endl << setprecision(5) << t << "\t\t" << 2.0*(pt-pl)/3.0 << endl;
	bulkplot << "t [fm]" << "\t\t" << "Pi [fm^-4]" << endl << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - B - p) << endl;
	plptplot << "t [fm]" << "\t\t" << "PL/PT" << endl << setprecision(5) << t << "\t\t" << (pl - B) / (pt - B) << endl;

	Tplot << "t [fm]" << "\t\t" << "T [fm^-1]" << endl << setprecision(5) << t << "\t\t" << T << endl;
	lambdaplot << "t [fm]" << "\t\t" << "lambda [fm^-1]" << endl << setprecision(5) << t << "\t\t" << lambda << endl;
	axplot << "t [fm]" << "\t\t" << "ax" << endl << setprecision(5) << t << "\t\t" << ax << endl;
	azplot << "t [fm]" << "\t\t" << "az" << endl << setprecision(5) << t << "\t\t" << az << endl;

	RpiInvplot << "t [fm]" << "\t\t" << "R_pi^-1" << endl << setprecision(5) << t << "\t\t" << sqrt(1.5) * 2.0*(pt-pl)/3.0 / p << endl;
	RbulkInvplot << "t [fm]" << "\t\t" << "R_Pi^-1" << endl << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - B - p) / p << endl;

	piNSplot << "t [fm]" << "\t\t" << "piNS" << endl << setprecision(5) << t << "\t\t" << piNS << endl;
	bulkNSplot << "t [fm]" << "\t\t" << "bulkNS" << endl << setprecision(5) << t << "\t\t" << bulkNS << endl;

	taupiplot << "t [fm]" << "\t\t" << "tau_pi" << endl << setprecision(5) << t << "\t\t" << taupi << endl;
	taubulkplot << "t [fm]" << "\t\t" << "tau_Pi" << endl << setprecision(5) << t << "\t\t" << taubulk << endl;

	Bplot << "t [fm]" << "\t\t" << "B" << endl << setprecision(5) << t << "\t\t" << B << endl;
	Beqplot << "t [fm]" << "\t\t" << "Beq" << endl << setprecision(5) << t << "\t\t" << Beq << endl;
	dBasyplot << "t [fm]" << "\t\t" << "dB2nd" << endl << setprecision(5) << t << "\t\t" << dBasy << endl;

	//evolution
	for(int i = 0; i < n; i++)
	{
		try{
			// compute intermediate values with Euler step
			Ttt_mid = Ttt + dt * dTtt_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Ttx_mid = Ttx + dt * dTtx_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Tty_mid = Tty + dt * dTty_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Ttn_mid = Ttn + dt * dTtn_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			pl_mid = pl + dt * dpl_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			pt_mid = pt + dt * dpt_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);
			B_mid = B + dt * dB_dt(Ttt,Ttx,Tty,Ttn,pl,pt,B,ut,ux,uy,un,e,p,lambda,ax,az,t);


			// find intermediate inferred and anisotropic variables
			get_inferred_variables(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,&ut,&ux,&uy,&un,&e,&p,t+dt);

			get_anisotropic_variables(e,pl_mid,pt_mid,B_mid,&lambda,&ax,&az);


			// add Euler step with respect to the intermediate value
			Ttt_end = Ttt_mid + dt * dTtt_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Ttx_end = Ttx_mid + dt * dTtx_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Tty_end = Tty_mid + dt * dTty_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Ttn_end = Ttn_mid + dt * dTty_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			pl_end = pl_mid + dt * dpl_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			pt_end = pt_mid + dt * dpt_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			B_end = B_mid + dt * dB_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,B_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);


			// increase time step
			t += dt;

			// Heun's Rule (average initial and end)
			// update variables at new time step
			Ttt = 0.5 * (Ttt + Ttt_end);
			Ttx = 0.5 * (Ttx + Ttx_end);
			Tty = 0.5 * (Tty + Tty_end);
			Ttn = 0.5 * (Ttn + Ttn_end);
			pl = 0.5 * (pl + pl_end);
			pt = 0.5 * (pt + pt_end);
			B = 0.5 * (B + B_end);


			// find inferred variables at new time step
			get_inferred_variables(Ttt,Ttx,Tty,Ttn,pl,pt,B,&ut,&ux,&uy,&un,&e,&p,t);
			get_anisotropic_variables(e,pl,pt,B,&lambda,&ax,&az);


			T = effectiveTemperature(e);

			Beq = equilibriumBquasi(T);

			m = T * z_Quasiparticle(T);

			piNS = 4.0 * (e+p) / (3.0*T*t) * shearViscosityToEntropyDensity(T);
			bulkNS = - (e+p) / (t*T) * bulkViscosityToEntropyDensity(T);

			// quasiparticle model
			taupi = (e+p) * shearViscosityToEntropyDensity(T) / (T*beta_shear(T));
			taubulk = (e+p) * bulkViscosityToEntropyDensity(T) / (T*beta_bulk(T));

			// dB2nd = -3.0*taubulk*mdmdT_Quasiparticle(T)/pow(z_Quasiparticle(T),2)*
			// 		speedOfSoundSquared(e)*(2.0*pt/3.0+pl/3.0-B-p)/(t*T);


			dBasy = -3.0*taubulk*mdmde_Quasiparticle(e)*(e+pl-B)*(2.0*pt/3.0+pl/3.0-B-p)/(t*m*m) /
			(1.0 + 4.0*taubulk*mdmde_Quasiparticle(e)*(e+pl-B)/(t*m*m));

			// write updated energy density to file
			if((i+1)%timesteps_per_write == 0)
			{
				// energy density, etc
				eplot << setprecision(5) << t << "\t\t" << e/e0 << "\t\t" << endl;
				piplot << setprecision(5) << t << "\t\t" << 2.0*(pt-pl)/3.0 << "\t\t" << endl;
				bulkplot << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - B - p) << "\t\t" << endl;
				plptplot << setprecision(5) << t << "\t\t" << (pl - B) / (pt - B) << "\t\t" << endl;

				// temperature and anisotropic variables
				Tplot << setprecision(5) << t << "\t\t" << T << "\t\t" << endl;
				lambdaplot << setprecision(5) << t << "\t\t" << lambda << "\t\t" << endl;
				axplot << setprecision(5) << t << "\t\t" << ax << "\t\t" << endl;
				azplot << setprecision(5) << t << "\t\t" << az << "\t\t" << endl;

				//cout << "Hey" << endl;

				//cout << "Lattice qcd energy:" << e << "		" << (I20_function(T,mbar) + equilibriumBquasi(T)) << endl;

				// inverse Reynolds numbers
				RpiInvplot << setprecision(5) << t << "\t\t" << sqrt(1.5) * 2.0*(pt-pl)/3.0 / p << "\t\t" << endl;
				RbulkInvplot << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - B - p) / p  << "\t\t" << endl;

				piNSplot << setprecision(5) << t << "\t\t" << piNS << "\t\t" << endl;
				bulkNSplot << setprecision(5) << t << "\t\t" << bulkNS  << "\t\t" << endl;


				taupiplot << setprecision(5) << t << "\t\t" << taupi << "\t\t" << endl;
				taubulkplot << setprecision(5) << t << "\t\t" << taubulk  << "\t\t" << endl;

				Bplot << setprecision(5) << t << "\t\t" << B  << "\t\t" << endl;
				Beqplot << setprecision(5) << t << "\t\t" << Beq  << "\t\t" << endl;
				dBasyplot << setprecision(5) << t << "\t\t" << dBasy  << "\t\t" << endl;
			}

		} catch (char const *excp)  {
        	cout << "\nRuntime error: " << excp;
        	exit(-1);
    	}
	}

	// close plot data files
	eplot.close();
	piplot.close();
	bulkplot.close();
	plptplot.close();

	Tplot.close();
	lambdaplot.close();
	axplot.close();
	azplot.close();

	RpiInvplot.close();
	RbulkInvplot.close();

	piNSplot.close();
	bulkNSplot.close();

	taupiplot.close();
	taubulkplot.close();

	Bplot.close();
	Beqplot.close();
	dBasyplot.close();



	printf("...done\n\n");


	duration = (clock() - begin) / (double)CLOCKS_PER_SEC;
    cout<<"Duration: " << setprecision(4) << duration / (double)n * 1000.0 << " ms per time step\n";


	return 0;
}






