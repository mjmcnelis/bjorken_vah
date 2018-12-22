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
const int gla_pts = 64;
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

 // 	int num_error;
	// // // Load gauss laguerre roots-weights
	// if((num_error = load_gauss_laguerre_data()) != 0)
	// {
	// 	fprintf(stderr, "Error loading gauss data (%d)!\n", num_error);
	// 	return 1;
	// }

	//for(int k = 0; k < gla_pts; k++) cout << setprecision(15) << root_gla[3][k] << ",";
	//printf("\n\n\n");
	//for(int k = 0; k < gla_pts; k++) cout << setprecision(15) << weight_gla[3][k] << ",";


    // int i;
    // for(i = 0; i < 10; i++)
    // {
    // 	if(i == 7)
    // 	{
    // 		cout << "Can see before break" << endl;
    // 		break;
    // 	}
    // 	if(i == 7)
    // 	{
    // 		cout << "Can't see after break" << endl;
    // 	}
    // }

    // cout << i << endl;

    //exit(-1);


	// input parameters
	double T0 = 0.5 * GEV_TO_INVERSE_FM;  // initial temperature in fm^-1
	double t0 = 0.25;				      // initial time in fm
	double tf = 100.0;					  // final time in fm




	// thermodynamic quantities
	double T = T0;
	double e0 = equilibriumEnergyDensity(T0);
	double e = e0;
	double p = equilibriumPressure(e);
	double p0 = p;
	double s = (e+p)/T;
	double cs2 = speedOfSoundSquared(e);
	double beq = equilibriumBquasi(T);
	double m = T * z_Quasiparticle(T);


	//cout << e << endl;
	//exit(-1);

	// viscosities and relaxation times
	double zetas = bulkViscosityToEntropyDensity(T);
	double etas = shearViscosityToEntropyDensity(T);
	double taupi = (s*etas) / beta_shear(T);
	double taubulk = (s*zetas) / beta_bulk(T);

	double piNS =  4.0 * (e0+p0) / (3.0*T0*t0) * shearViscosityToEntropyDensity(T0);
	double bulkNS = - (e0+p0) / (t0*T0) * bulkViscosityToEntropyDensity(T0);


	// initial Tmunu components
	double Ttt = e;
	double Ttx = 0.0;
	double Tty = 0.0;
	double Ttn = 0.0;

	// initial flow profile
	double ut = 1.0;
	double ux = 0.0;
	double uy = 0.0;
	double un = 0.0;

	// Glasma initial conditions
	double plptratio = 0.02;
	double pt = (3.0/(2.0+plptratio)) * e / 3.0;
	double pl = (3.0 - 6.0/(2.0+plptratio)) * e / 3.0;


	double dbasy = -3.0*taubulk*mdmde_Quasiparticle(e)*(e+pl)*(2.0*pt/3.0+pl/3.0-p)/(t0*m*m) /
			(1.0 + 4.0*taubulk*mdmde_Quasiparticle(e)*(e+pl)/(t0*m*m));

	//cout << beq + dbasy << endl;

	double b_default = (beq + dbasy);

	//cout << "pl_kin = " << pl + b << endl;
	//cout << "pt_kin = " << pt + b << endl;
	//cout << "e_kin = " << e - b << endl;

	double lambda = T;
	double ax = 1.0;
	double az = 1.0;
	double b;
	if(plptratio > 0.08) b = b_default;
	else b = b_default * (0.14 + (1.0 - 0.14) * (plptratio - 0.01) / 0.07);

	//b = 0.14 * b_default;

	//cout << 0.2 * beq / b_default << endl;

	// int N = 100;


	// double frac_first = 1.0;
	// double frac_last = 1.0;

	// int default_counter = 1;	// assume default b solves equation

	// // try solving default equation
	// try
	// {
	// get_anisotropic_variables(e, pl, pt, b_default, &lambda, &ax, &az);
	// }
	// catch (char const *excp)
	// {
 //        	cout << "\nInitialization error: " << excp;
 //        	cout << endl << "Increasing mean field..." << endl;
 //        	default_counter = 0;
 //    }

 //    // increase b until 1st solution found
 //    if(!default_counter)
 //    {
 //    	int modb_counter = 1;	// assume modified b solves equation
 //    	//int modb_min_counter = 0;
 //    	int modb_1st_counter = 0;

 //    	//double lambda = T;
 //    	//double ax = 1.0;
 //    	//double az = 1.0;

 //    	for(int i = 1; i < N; i++)
	// 	{
	// 		double del = 1.0 / (double) N;
	// 		double frac = 1.0 - (double)i * del;

	// 		double b = frac * b_default;

	// 		double lambda = T;
	// 		double ax = 1.0;
	// 		double az = 1.0;

	// 		// I can make it faster if I kept track of the updated variables
	// 		// why did it crash?


	// 		// note: I want the two fracs to be different by a fair amount so I might want to decrease N. Otherwise the code might not run.

	// 		try
	// 		{
	// 			modb_counter = 1; // try assumed solution
	// 			get_anisotropic_variables(e, pl, pt, b, &lambda, &ax, &az);


	// 		}
	// 		catch (char const *excp)
	// 		{
	// 			//cout << "\nInitialization error: " << excp;
	// 	        modb_counter = 0; // not the solution
	// 	    }

	// 	    if(modb_counter && modb_1st_counter == 0)
	// 	    {
	// 	    	modb_1st_counter = 1;
	// 	    	cout << "Found 1st solution at frac = " << setprecision(5) << frac << endl;
	// 	    	frac_first = frac - del;
	// 	    }

	// 	    if(modb_counter == 0 && modb_1st_counter == 1)
	// 	    {
	// 	    	cout << "Found last solution at frac = " << setprecision(5) << frac + del << endl;

	// 	    	frac_last = frac + del;

	// 	    	break;
	// 	    }
	// 	}
 //    }

	// lambda = T;
	// ax = 1.0;
	// az = 1.0;

	// //double b = frac_first * b_default;
	// double b = frac_last * b_default;

	try
	{
	get_anisotropic_variables(e, pl, pt, b, &lambda, &ax, &az);
	}
	catch (char const *excp)
	{
        	cout << "\nInitialization error: " << excp;
        	exit(-1);
    }

    cout << "\nT = " << T << endl;
	cout << "lambda = " << setprecision(4) << lambda / GEV_TO_INVERSE_FM << " * GEV_TO_INVERSE_FM;" << endl;
	cout << "ax = " << ax << ";" << endl;
	cout << "az = " << az << ";" <<endl;
	cout << "b/beq = " << b / beq << endl;

	//exit(-1);


	// intermediate and end values
	double Ttt_mid, Ttt_end;
	double Ttx_mid, Ttx_end;
	double Tty_mid, Tty_end;
	double Ttn_mid, Ttn_end;
	double pl_mid, pl_end;
	double pt_mid, pt_end;
	double b_mid, b_end;


	// time data (t = tau)
	double t = t0;
	const double dt = 0.001;
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
	bulkplot << "t [fm]" << "\t\t" << "Pi [fm^-4]" << endl << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - p) << endl;
	plptplot << "t [fm]" << "\t\t" << "PL/PT" << endl << setprecision(5) << t << "\t\t" << pl/pt << endl;

	Tplot << "t [fm]" << "\t\t" << "T [fm^-1]" << endl << setprecision(5) << t << "\t\t" << T << endl;
	lambdaplot << "t [fm]" << "\t\t" << "lambda [fm^-1]" << endl << setprecision(5) << t << "\t\t" << lambda << endl;
	axplot << "t [fm]" << "\t\t" << "ax" << endl << setprecision(5) << t << "\t\t" << ax << endl;
	azplot << "t [fm]" << "\t\t" << "az" << endl << setprecision(5) << t << "\t\t" << az << endl;

	RpiInvplot << "t [fm]" << "\t\t" << "R_pi^-1" << endl << setprecision(5) << t << "\t\t" << sqrt(1.5) * 2.0*(pt-pl)/3.0 / p << endl;
	RbulkInvplot << "t [fm]" << "\t\t" << "R_Pi^-1" << endl << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - p) / p << endl;

	piNSplot << "t [fm]" << "\t\t" << "R_piNS^-1" << endl << setprecision(5) << t << "\t\t" << sqrt(1.5)*piNS/p << endl;
	bulkNSplot << "t [fm]" << "\t\t" << "R_bulkNS^-1" << endl << setprecision(5) << t << "\t\t" << bulkNS/p << endl;

	taupiplot << "t [fm]" << "\t\t" << "tau_pi" << endl << setprecision(5) << t << "\t\t" << taupi << endl;
	taubulkplot << "t [fm]" << "\t\t" << "tau_Pi" << endl << setprecision(5) << t << "\t\t" << taubulk << endl;

	Bplot << "t [fm]" << "\t\t" << "B" << endl << setprecision(5) << t << "\t\t" << b << endl;
	Beqplot << "t [fm]" << "\t\t" << "Beq" << endl << setprecision(5) << t << "\t\t" << beq << endl;
	dBasyplot << "t [fm]" << "\t\t" << "dB2nd" << endl << setprecision(5) << t << "\t\t" << dbasy << endl;

	//evolution
	for(int i = 0; i < n; i++)
	{
		try{
			// compute intermediate values with Euler step
			Ttt_mid = Ttt + dt * dTtt_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Ttx_mid = Ttx + dt * dTtx_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Tty_mid = Tty + dt * dTty_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			Ttn_mid = Ttn + dt * dTtn_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			pl_mid = pl + dt * dpl_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			pt_mid = pt + dt * dpt_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);
			b_mid = b + dt * db_dt(Ttt,Ttx,Tty,Ttn,pl,pt,b,ut,ux,uy,un,e,p,lambda,ax,az,t);


			// find intermediate inferred and anisotropic variables
			get_inferred_variables(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,&ut,&ux,&uy,&un,&e,&p,t+dt);

			get_anisotropic_variables(e,pl_mid,pt_mid,b_mid,&lambda,&ax,&az);


			// add Euler step with respect to the intermediate value
			Ttt_end = Ttt_mid + dt * dTtt_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Ttx_end = Ttx_mid + dt * dTtx_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Tty_end = Tty_mid + dt * dTty_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			Ttn_end = Ttn_mid + dt * dTty_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			pl_end = pl_mid + dt * dpl_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			pt_end = pt_mid + dt * dpt_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);
			b_end = b_mid + dt * db_dt(Ttt_mid,Ttx_mid,Tty_mid,Ttn_mid,pl_mid,pt_mid,b_mid,ut,ux,uy,un,e,p,lambda,ax,az,t+dt);


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
			b = 0.5 * (b + b_end);


			// find inferred variables at new time step
			get_inferred_variables(Ttt,Ttx,Tty,Ttn,pl,pt,b,&ut,&ux,&uy,&un,&e,&p,t);
			get_anisotropic_variables(e,pl,pt,b,&lambda,&ax,&az);


			T = effectiveTemperature(e);
			cs2 = speedOfSoundSquared(e);
			s = (e+p)/T;
			beq = equilibriumBquasi(T);
			m = T * z_Quasiparticle(T);

			piNS = 4.0 * (e+p) / (3.0*T*t) * shearViscosityToEntropyDensity(T);
			bulkNS = - (e+p) / (t*T) * bulkViscosityToEntropyDensity(T);

			// quasiparticle model
			taupi = s * shearViscosityToEntropyDensity(T) / beta_shear(T);
			taubulk = s * bulkViscosityToEntropyDensity(T) / beta_bulk(T);

			//dBasy = -3.0*taubulk*(cs2/s)*mdmdT_Quasiparticle(T)*(e+pl-B)*(2.0*pt/3.0+pl/3.0-B-p)/(t*m*m) /
			//(1.0 + 4.0*taubulk*(cs2/s)*mdmdT_Quasiparticle(T)*(e+pl-B)/(t*m*m));

			// dBasy = -3.0*taubulk*mdmde_Quasiparticle(e)*(e+pl-B)*(2.0*pt/3.0+pl/3.0-B-p)/(t*m*m) /
			// (1.0 + 4.0*taubulk*mdmde_Quasiparticle(e)*(e+pl-B)/(t*m*m));

			dbasy = -3.0*taubulk*mdmde_Quasiparticle(e)*(e+pl)*(2.0*pt/3.0+pl/3.0-p)/(t*m*m) /
			(1.0 + 4.0*taubulk*mdmde_Quasiparticle(e)*(e+pl)/(t*m*m));

			// write data to file
			if((i+1)%timesteps_per_write == 0)
			{
				// energy density, etc
				eplot << setprecision(5) << t << "\t\t" << e/e0 << "\t\t" << endl;
				piplot << setprecision(5) << t << "\t\t" << 2.0*(pt-pl)/3.0 << "\t\t" << endl;
				bulkplot << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - p) << "\t\t" << endl;
				plptplot << setprecision(5) << t << "\t\t" << pl/pt << "\t\t" << endl;

				// temperature and anisotropic variables
				Tplot << setprecision(5) << t << "\t\t" << T << "\t\t" << endl;
				lambdaplot << setprecision(5) << t << "\t\t" << lambda << "\t\t" << endl;
				axplot << setprecision(5) << t << "\t\t" << ax << "\t\t" << endl;
				azplot << setprecision(5) << t << "\t\t" << az << "\t\t" << endl;

				// inverse Reynolds numbers
				RpiInvplot << setprecision(5) << t << "\t\t" << sqrt(1.5) * 2.0*(pt-pl)/3.0 / p << "\t\t" << endl;
				RbulkInvplot << setprecision(5) << t << "\t\t" << (2.0*pt/3.0 + pl/3.0 - p) / p  << "\t\t" << endl;

				piNSplot << setprecision(5) << t << "\t\t" << sqrt(1.5)*piNS/p << "\t\t" << endl;
				bulkNSplot << setprecision(5) << t << "\t\t" << bulkNS/p  << "\t\t" << endl;

				taupiplot << setprecision(5) << t << "\t\t" << taupi << "\t\t" << endl;
				taubulkplot << setprecision(5) << t << "\t\t" << taubulk  << "\t\t" << endl;

				Bplot << setprecision(5) << t << "\t\t" << b  << "\t\t" << endl;
				Beqplot << setprecision(5) << t << "\t\t" << beq  << "\t\t" << endl;
				dBasyplot << setprecision(5) << t << "\t\t" << dbasy  << "\t\t" << endl;
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

	printf("\n...done\n\n");

	duration = (clock() - begin) / (double)CLOCKS_PER_SEC;
    cout<<"Duration: " << setprecision(4) << duration / (double)n * 1000.0 << " ms per time step\n";


	return 0;
}






