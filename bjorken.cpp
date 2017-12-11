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
	const double T0 = 0.6 * GEV_TO_INVERSE_FM;  // initial temperature in fm^-1
	const double tau0 = 0.25;					// initial time in fm
	const double tauf = 30.0;					// final time in fm

	double mbar0 = z_Quasiparticle(T0);

	// initial anisotropic parameters (equilibrium initial conditions)
	double lambda = T0;
	double ax = 1.0;
	double az = 1.0;

	double T = T0;    // for plotting purposes
	double mbar = mbar0;
	//double pkinetic = equilibriumKineticPressure(T0);



	double pkinetic = I21_function(T0,mbar0);



	// initial flow profile
	double ut = 1.0;
	double ux = 0.0;
	double uy = 0.0;
	double un = 0.0;


	// initial energy density and pressure  (units = [fm^-4])
	const double e0 = equilibriumEnergyDensity(T0);
	double e = e0;
	const double p0 = equilibriumPressure(e0);
	double p = p0;

	double piNS =  4.0 * (e0+p0) / (3.0*T0*tau0) * shearViscosityToEntropyDensity(T0);
	double bulkNS = - (e0+p0) / (tau0*T0) * bulkViscosityToEntropyDensity(T0);



	// initial T^{\tau\mu} components (units = [fm^-4])
	double Ttt = e0;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttn = 0.0;

	// initialize kinetic pl and pt (quasiparticle equilibrium w/o B(T))
	//double pl = equilibriumKineticPressure(T0);
	//double pt = equilibriumKineticPressure(T0);   // temporary

	double pl = I21_function(T0,mbar0);
	double pt = I21_function(T0,mbar0);


	// intermediate and end values
	double Ttt_mid, Ttt_end;
	double Ttx_mid, Ttx_end;
	double Tty_mid, Tty_end;
	double Ttn_mid, Ttn_end;
	double pl_mid, pl_end;
	double pt_mid, pt_end;


	// time data
	double tau = tau0;
	const double dtau = 0.001;
	const int n = floor((tauf - tau0) / dtau);
	const int timesteps_per_write = 10;

	// Data files for plots
	ofstream eplot, piplot, bulkplot, plptplot;
	ofstream Tplot, lambdaplot, axplot, azplot;
	ofstream RpiInvplot, RbulkInvplot;
	ofstream piNSplot, bulkNSplot;

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

	eplot << "tau [fm]" << "\t\t" << "e/e0" << endl << setprecision(5) << tau << "\t\t" << e/e0 << endl;
	piplot << "tau [fm]" << "\t\t" << "pi [fm^-4]" << endl << setprecision(5) << tau << "\t\t" << 0.0 << endl;
	bulkplot << "tau [fm]" << "\t\t" << "Pi [fm^-4]" << endl << setprecision(5) << tau << "\t\t" << 0.0 << endl;
	plptplot << "tau [fm]" << "\t\t" << "PL/PT" << endl << setprecision(5) << tau << "\t\t" << (p + pl - equilibriumKineticPressure(T)) / (p + pt - equilibriumKineticPressure(T)) << endl;

	Tplot << "tau [fm]" << "\t\t" << "T [fm^-1]" << endl << setprecision(5) << tau << "\t\t" << T << endl;
	lambdaplot << "tau [fm]" << "\t\t" << "lambda [fm^-1]" << endl << setprecision(5) << tau << "\t\t" << lambda << endl;
	axplot << "tau [fm]" << "\t\t" << "ax" << endl << setprecision(5) << tau << "\t\t" << ax << endl;
	azplot << "tau [fm]" << "\t\t" << "az" << endl << setprecision(5) << tau << "\t\t" << az << endl;

	RpiInvplot << "tau [fm]" << "\t\t" << "R_pi^-1" << endl << setprecision(5) << tau << "\t\t" << 0.0 << endl;
	RbulkInvplot << "tau [fm]" << "\t\t" << "R_Pi^-1" << endl << setprecision(5) << tau << "\t\t" << 0.0 << endl;

	piNSplot << "tau [fm]" << "\t\t" << "piNS" << endl << setprecision(5) << tau << "\t\t" << piNS << endl;
	bulkNSplot << "tau [fm]" << "\t\t" << "bulkNS" << endl << setprecision(5) << tau << "\t\t" << bulkNS << endl;


	// evolution
	for(int i = 0; i < n; i++)
	{
		// compute intermediate values with Euler step
		Ttt_mid = Ttt + dtau * dTtt_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);
		Ttx_mid = Ttx + dtau * dTtx_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);
		Tty_mid = Tty + dtau * dTty_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);
		Ttn_mid = Ttn + dtau * dTtn_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);
		pl_mid = pl + dtau * dpl_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);
		pt_mid = pt + dtau * dpt_dtau(Ttt, Ttx, Tty, Ttn, pl, pt, ut, ux, uy, un, e, p, lambda, ax, az, tau);




		// find intermediate inferred and anisotropic variables
		get_inferred_variables(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, &ut, &ux, &uy, &un, &e, &p, tau + dtau);

		get_anisotropic_variables(e, pl_mid, pt_mid, &lambda, &ax, &az);



		// add Euler step with respect to the intermediate value
		Ttt_end = Ttt_mid + dtau * dTtt_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);
		Ttx_end = Ttx_mid + dtau * dTtx_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);
		Tty_end = Tty_mid + dtau * dTty_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);
		Ttn_end = Ttn_mid + dtau * dTty_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);
		pl_end = pl_mid + dtau * dpl_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);
		pt_end = pt_mid + dtau * dpt_dtau(Ttt_mid, Ttx_mid, Tty_mid, Ttn_mid, pl_mid, pt_mid, ut, ux, uy, un, e, p, lambda, ax, az, tau + dtau);

		// increase time step
		tau += dtau;


		// Heun's Rule (average initial and end)
		// update variables at new time step
		Ttt = 0.5 * (Ttt + Ttt_end);
		Ttx = 0.5 * (Ttx + Ttx_end);
		Tty = 0.5 * (Tty + Tty_end);
		Ttn = 0.5 * (Ttn + Ttn_end);
		pl = 0.5 * (pl + pl_end);
		pt = 0.5 * (pt + pt_end);


		// find inferred variables at new time step
		get_inferred_variables(Ttt, Ttx, Tty, Ttn, pl, pt, &ut, &ux, &uy, &un, &e, &p, tau);

		get_anisotropic_variables(e, pl, pt, &lambda, &ax, &az);

		T = effectiveTemperature(e);  // for plots
		//pkinetic = equilibriumKineticPressure(T);

		mbar = z_Quasiparticle(T);
		pkinetic = I21_function(T,mbar);

		piNS = 4.0 * (e+p) / (3.0*T*tau) * shearViscosityToEntropyDensity(T);
		bulkNS = - (e+p) / (tau*T) * bulkViscosityToEntropyDensity(T);


		// write updated energy density to file
		if((i+1)%timesteps_per_write == 0)
		{
			// energy density, etc
			eplot << setprecision(5) << tau << "\t\t" << e/e0 << "\t\t" << endl;
			piplot << setprecision(5) << tau << "\t\t" << 2.0*(pt-pl)/3.0 << "\t\t" << endl;
			bulkplot << setprecision(5) << tau << "\t\t" << (2.0*pt/3.0 + pl/3.0 - pkinetic) << "\t\t" << endl;
			plptplot << setprecision(5) << tau << "\t\t" << (p + pl - pkinetic) / (p + pt - pkinetic) << "\t\t" << endl;

			// temperature and anisotropic variables
			Tplot << setprecision(5) << tau << "\t\t" << T << "\t\t" << endl;
			lambdaplot << setprecision(5) << tau << "\t\t" << lambda << "\t\t" << endl;
			axplot << setprecision(5) << tau << "\t\t" << ax << "\t\t" << endl;
			azplot << setprecision(5) << tau << "\t\t" << az << "\t\t" << endl;

			//cout << "Hey" << endl;

			// inverse Reynolds numbers
			RpiInvplot << setprecision(5) << tau << "\t\t" << sqrt(1.5) * 2.0*(pt-pl)/3.0 / p << "\t\t" << endl;
			RbulkInvplot << setprecision(5) << tau << "\t\t" << (2.0*pt/3.0 + pl/3.0 - pkinetic) / p  << "\t\t" << endl;

			piNSplot << setprecision(5) << tau << "\t\t" << piNS << "\t\t" << endl;
			bulkNSplot << setprecision(5) << tau << "\t\t" << bulkNS  << "\t\t" << endl;
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

	printf("...done\n\n");


	duration = (clock() - begin) / (double)CLOCKS_PER_SEC;
    cout<<"Duration: " << setprecision(4) << duration / (double)n * 1000.0 << " ms per time step\n";


	return 0;
}






