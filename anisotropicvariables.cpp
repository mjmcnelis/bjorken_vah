
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <exception>
using namespace std;
#include "anisotropicvariables.hpp"
#include "qcd.hpp"
#include "anisotropic_integrands.hpp"
#include "gauss_integration.hpp"
#include "anisotropic_transport.hpp"
#define EPS_MIN 1.0e-16


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                    LUP LINEAR EQUATION SOLVER                    ::
//                                                                  ::
//     Solves linear equation Ax = b using LU decomposition         ::
//	   with implicit partial pivoting (LUP). To directly solve      ::
//     for Ax = b, run these two functions concurrently:            ::
//                                                                  ::
//            LUP_decomposition              LUP_solve              ::
//                                                                  ::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void LUP_decomposition(double ** A, int n, int * pvector)
{
	// A = n x n matrix; function does A -> PA = LU; (L,U) of PA stored in same array
	// n = size of A
	// pvector = permutation vector; set initial pvector[i] = i (to track implicit partial pivoting)
	// 			 function updates pvector if there are any rows exchanges in A (to be used on b in LUP_solve)

	int i;	   // rows
	int j;	   // columns
	int k;     // dummy matrix index
	int imax;  // pivot row index
	double big;
	double sum;
	double temp;
	double implicit_scale[n];

	// Initialize permutation vector
	// to default no-pivot values
	for(i = 0; i < n; i++)
	{
		pvector[i] = i;
	}
	// Implicit scaling info. for A
	for(i = 0; i < n; i++)
	{
		big = 0.0;
		for(j = 0; j < n; j++)
		{
			temp = fabs(A[i][j]);
			if(temp > big)
			{
				big = temp;  // update biggest element in the ith row
			}
		}
		if(big == 0.0)
		{
			printf("Singular matrix in the routine");
			break;
		}
		implicit_scale[i] = 1.0 / big;  // store implicit scale of row i (will be used later)
	}
	// LU Decomposition
	for(j = 0; j < n; j++)
	{
		// loop rows i = 1 to j - 1
		for(i = 0; i < j; i++)
		{
			sum = A[i][j];
			for(k = 0; k < i; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;  // update U[i][j] elements
		}

		big = 0.0;          // initialize search for the largest normalized pivot in j column

		// loop through rows i = j to n-1
		for(i = j; i < n; i++)
		{
			sum = A[i][j];
			for(k = 0; k < j; k++)
			{
				sum -= A[i][k] * A[k][j];
			}
			A[i][j] = sum;   // update U[j][j] and L[i][j] elements (no division in L yet until the pivot determined)

			temp = implicit_scale[i] * fabs(sum);
			if(temp >= big)
			{
				big = temp;  // searchs for the biggest normalized member (imax) in column j
				imax = i;	 // implicit scale * A[i][j] normalizes each row entry i in column j before comparing
			}			     // implicit scale different for each i, that's why it's important
		}
		if(j != imax)
		{
			// then exchange rows j and imax of A
			for(k = 0; k < n; k++)
				{
					temp = A[imax][k];
					A[imax][k] = A[j][k];
					A[j][k] = temp;
				}
			implicit_scale[imax] = implicit_scale[j];   // interchange scaling
		}

		pvector[j] = imax;   		  // update permutation vector keeps track of pivot indices

		if(A[j][j] == 0.0)
		{
			A[j][j] = EPS_MIN;        // matrix is singular
		}
		if(j != n-1)                  // there is no L[n,n] element
		{
			temp = 1.0 / A[j][j];     // divide L[i,j] elements by the pivot
			for(i = j+1; i < n; i++)
			{
				A[i][j] *= temp;
			}
		}
	}
}

void LUP_solve(double ** PA, int n, int * pvector, double b[])
{
	// input vector b is transformed to the solution x  of Ax = b
	// PA is the input permutated matrix from LUP_decomposition (will not be updated here)
	// input pvector comes from LUP_decomposition (generally not default); used to switch rows of b
	int i;       // rows
	int j;       // columns
	int m = -1;  // used to skip first few for loops (j) involving a zero (priorly permutated) b[j] element (pointless to multiply by 0)
	int ip;      // assigned permutation row index pvector[i]
	double sum;
	// Forward substitution routine for Ly = b
	for(i = 0; i < n; i++)
	{
		ip = pvector[i];         // permute b accordingly
		sum = b[ip];             // starting value given right b[ip]
		b[ip] = b[i];            // switch value of b[ip] to b[i]
		if(m != -1)			     // if m stays -1, skip for loop until it's been assigned a j_min >= 0
		{
			for(j = m; j <= i-1; j++)
			{
				sum -= PA[i][j] * b[j];    // forward iteration
			}
		}
		else if(sum)
		{
			m = i;               // once encounter nonzero b[ip], m = j_min >= 0 from now on
		}
		b[i] = sum;              // update y[i] and store in b
	}
	// Backward substitution routine for Ux = y
	for(i = n-1; i >= 0; i--)
	{
		sum = b[i];
		for(j = i+1; j < n; j++)
		{
			sum -= PA[i][j] * b[j];       // backward iteration
		}
		b[i] = sum / PA[i][i];            // update x[i] and store in b (FINAL SOLUTION)
	}
}

void free_2D(double ** M, int n)
{
	for (int i = 0; i < n; i++) free(M[i]);
    free(M);
}












//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                       ANISOTROPIC VARIABLES                      ::
//                                                                  ::
//     Compute the anisotropic variables (lambda,ax,az) from        ::
//	   the quantities (e,B,pl,pt)       		                    ::
//                                                                  ::
//                     get_anisotropic_variables                    ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::






void get_anisotropic_variables(double e, double pl, double pt, double B, double *lambda, double *ax, double *az)
{
	// not designed for conformal mode, it will crash

	// (e,B,pt,pl) should already be updated before running this
	// order: update conserved variables, update inferred variables, update anisotropic variables

	const double T = effectiveTemperature(e);				   // temperature


	//const double Ekinetic = equilibriumKineticEnergyDensity(T);    // quasiparticle kinetic energy density

	//const double mbareq = z_Quasiparticle(T);
	//const double Ekinetic = I20_function(T,mbareq);
	//const double dB = B - equilibriumBquasi(T);

	//const double Ea = Ekinetic - dB;
	const double Ea = e - B;
	const double PTa = pt;
	const double PLa = pl;

	//throw "Boom ";

	if(Ea < 0.0)
	{
		throw "Ea is out of bounds!\n";
	}
	if(PTa < 0.0)
	{
		throw "PTa is out of bounds!\n";
	}
	if(PLa < 0.0)
	{
		throw "PLa is out of bounds!\n";
	}

	const double g = 51.4103536012791;                       // degeneracy factor g (nf = 3 flavors)

	const int pbar_pts = 32;

	double pbar_rootF[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

	double pbar_weightF[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

	double pbar_rootJ[pbar_pts] = {0.299618729049241,0.701981065353977,1.24974569814569,1.94514382443706,2.78994869155499,3.78614305529879,4.93605293430173,6.24240692441721,7.70838362900271,9.3376617765878,11.1344784990019,13.1036993239838,15.2509033992287,17.5824881879873,20.1057991514724,22.8292918399089,25.7627366009885,28.9174802233293,32.3067850039226,35.9462752181738,39.8545359681502,44.0539338428991,48.5717701879893,53.4419507581545,58.7074908793654,64.4244418290391,70.6683893377061,77.5459900633602,85.2174664086695,93.9467116599065,104.238552969691,117.391742318923};

	double pbar_weightJ[pbar_pts] = {0.00660146448073508,0.0813584931042281,0.347537436309438,0.809963198105261,1.22739584119905,1.32050782861975,1.06049919505728,0.655616488144915,0.318173017008472,0.122743109012855,0.0379333897858022,0.00943187028689987,0.00188978713293874,0.000304914974586437,3.95130877631855e-05,4.09377958251348e-06,3.36921618654073e-07,2.1841295448875e-08,1.10337736506627e-09,4.28638379146177e-11,1.25966453444067e-12,2.74423030367617e-14,4.32175197361363e-16,4.76686817705967e-18,3.53643350342934e-20,1.67355018349782e-22,4.70254099995936e-25,7.09116556196869e-28,4.93082516196282e-31,1.23284946609868e-34,6.91389702736573e-39,2.63586492716958e-44};


	double lambdai = *lambda;			                     // initial guess for anisotropic variables
	double axi = *ax;										 // are previous time step values
	double azi = *az;


	const int n = 3; 										 // dimension space
	double X[n] = {lambdai, axi, azi};						 // initial guess vector
	double Xold[n] = {lambdai,axi,azi};						 // original guess
  	double F[n];  											 // F vector (root equation: F[X] = 0)
	// J = Jacobian of F
	double **J = (double **) malloc(n * sizeof(double *));
	for(int k = 0; k < n; k++) J[k] = (double *) malloc(n* sizeof(double));
 	int pvector[n];									  		 // permutation vector


 	// mbar = m(T) / lambda
	double thermal_mass = z_Quasiparticle(T) * T;	     // m(T) fixed wpt lambda
	double mbari;                   				     // mbar = m(T) / lambda


	// prefactors and factors of F/J elements: to be evaluated at ith iteration of X
 	double commonfactori;
  	double lambdai2, lambdai3, axi2, azi2, lambdaiaxi3, lambdaiazi3;
  	double factorEai, factorPTai, factorPLai;
  	double factorI2001, factorI2011, factorI2201, factorI402m1, factorI421m1, factorI440m1;

  	double factorI401m1, factorI420m1;


	// anisotropic functions evaluated at ith iteration of X
	double Eai, PTai, PLai;
	double I2001, I2011, I2201, I402m1, I421m1, I440m1;

	double I401m1, I420m1;



	int i = 0;				  // starting ith iteration
	int Nmax = 1000;	      // max number of iterations
	double dXnorm2;           // L2-norm of dX iteration
	double Fnorm2;		      // L2-norm of F
	double toldX = 1.0e-8;    // tolerance for dX
	double tolF = 1.0e-11;    // tolerance for F

	toldX = 1.0e-7;
	tolF = 1.0e-10;

	double l = 0.1;

	int m = 0;
	// I should spend Christmas break working on the Broydyn method

	// Find anisotropic variables using 3D Newton Method
	do{

		//cout << i << ", ";

		// Evaluate mass parameter
		mbari = thermal_mass / lambdai;

		// Evaluate factors and prefactors
		lambdai2 = lambdai * lambdai;
	    lambdai3 = lambdai2 * lambdai;
	    axi2 = axi * axi;
	    azi2 = azi * azi;
	    lambdaiaxi3 = lambdai * axi * axi2;
	    lambdaiazi3 = lambdai * azi * azi2;

	    // Default prefactor for the momemt (n,l,q,s) = 0
	    commonfactori = g * axi2 * azi * lambdai2 / (4.0*M_PI*M_PI);

	    factorEai = commonfactori * lambdai2;
	    factorPTai = commonfactori * axi2 * lambdai2 / 2.0;
	    factorPLai = commonfactori * azi2 * lambdai2;

	    factorI2001 = commonfactori * lambdai3;
	    factorI2011 = commonfactori * axi2 * lambdai3 / 2.0;
	    factorI2201 = commonfactori * azi2 * lambdai3;

	    factorI401m1 = factorI2011;
	    factorI420m1 = factorI2201;

	    factorI402m1 = commonfactori * axi2 * axi2 * lambdai3 / 8.0;
	    factorI421m1 = commonfactori * axi2 * azi2 * lambdai3 / 2.0;
	    factorI440m1 = commonfactori * azi2 * azi2 * lambdai3;



	    // Evaluate anisotropic functions for F  (1D Gauss Laguerre integrals)
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	    Eai = factorEai * Gauss_Aniso_1D(Ea_integrand, pbar_rootF, pbar_weightF, pbar_pts, axi, azi, mbari);
	    PTai = factorPTai * Gauss_Aniso_1D(PTa_integrand, pbar_rootF, pbar_weightF, pbar_pts, axi, azi, mbari);
	    PLai = factorPLai * Gauss_Aniso_1D(PLa_integrand, pbar_rootF, pbar_weightF, pbar_pts, axi, azi, mbari);


	    ////////////////////////////
	    //                        //
	    //   F  =  Eai - Ea       //
	    //         PTai - PTa     //
	    //         PLai - PLa     //
	    //                        //
	    ////////////////////////////


		F[0] = Eai - Ea;
    	F[1] = PTai - PTa;
    	F[2] = PLai - PLa;

		// Calculate L2-norm of F
    	Fnorm2 = sqrt(fabs(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]));

    	// Evaluate anisotropic functions for J (more 1D Gauss-Laguerre integrals)
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	    I2001 = factorI2001 * Gauss_Aniso_1D(I2001_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);
	    I2011 = factorI2011 * Gauss_Aniso_1D(I2011_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);
	    I2201 = factorI2201 * Gauss_Aniso_1D(I2201_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);

	    I401m1 = factorI401m1 * Gauss_Aniso_1D(I401m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);
	    I420m1 = factorI420m1 * Gauss_Aniso_1D(I420m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);

	    I402m1 = factorI402m1 * Gauss_Aniso_1D(I402m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);
	    I421m1 = factorI421m1 * Gauss_Aniso_1D(I421m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);
	    I440m1 = factorI440m1 * Gauss_Aniso_1D(I440m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, axi, azi, mbari);


	    //////////////////////////////////////////////////////////////////////////////
	    //                                                                          //
	    //    J  =  I2001/lambda2   2*I401m1/lambda/ax3    I420m1/lambda/az3        //
	    //                                                                          //
	    //          I2011/lambda2   4*I402m1/lambda/ax3    I421m1/lambda/az3        //
	    //                                                                          //
	    //          I2201/lambda2   2*I421m1/lambda/ax3    I440m1/lambda/az3        //
	    //                                                                          //
	    //////////////////////////////////////////////////////////////////////////////

	    // row 1
	    J[0][0] = I2001/lambdai2;
	    J[0][1] = 2.0 * I401m1 / lambdaiaxi3;
	    //J[0][1] = 2.0*(Eai+PTai)/axi;
	    J[0][2] = I420m1 / lambdaiazi3;
	    //J[0][2] = (Eai+PLai)/azi;
	    // row 2
	    J[1][0] = I2011/lambdai2;
	    J[1][1] = 4.0 * I402m1 / lambdaiaxi3;
	    J[1][2] = I421m1 / lambdaiazi3;
	    // row 3
	    J[2][0] = I2201/lambdai2;
	    J[2][1] = 2.0 * I421m1 / lambdaiaxi3;
	    J[2][2] = I440m1 / lambdaiazi3;



	    // Solve matrix equation: J * dX = - F
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	    for(int k = 0; k < n; k++) F[k] = - F[k];  // change sign of F first
	    // LU solver routine
	    LUP_decomposition(J, n, pvector);          // LUP of J now stored in J (pvector also updated)
	    LUP_solve(J, n, pvector, F);               // dX is now stored in F
	    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	    // here F is the full Newton step dX

	   //  l = 1.0;

	   //  // check whether or not F (or f?) decreased:
	   //  // it'll go something like
	   //  if(f > fold)
	   //  {
	   //  	// if 2nd:Nmax Newton fails use cubic g(l) solution
	   //  	if(i > 0)
	   //  	{
	   //  		// work out all details of the cubic solution

	   //  		g0 = F0;
	   //  		g1 = F1;     // these guys need some work...
	   //  		g2 = F2;
	   //  		gprime0 = gradf * dX;

	   //  		a = ((g1-gprime0*l1-g0)/(l1*l1) - (g2-gprime0*l2-g0)/(l2*l2)) / (l1-l2);
				// b = (-l2*(g1-gprime0*l1-g0)/(l1*l1) + l1*(g2-gprime0*l2-g0)/(l2*l2)) / (l1-l2);

				// // update current, previous and 2nd previous l's
				// l2 = l1;
				// l1 = l;
				// l = (-b + sqrt(b*b - 3.0*a*gprime0)) / (3.0*a);

				// if(l < 0.1*l2) // l1 used in routine stored in l2
				// 	l = 0.1 * l2;
				// if(l > 0.5*l2)
				// 	l = 0.5*l2;
	   //  	}
	   //  	// if 1st Newton fails use quadratic g(l) solution
	   //  	else if(i == 0)
	   //  	{
	   //  		// work out details of the quadratic solution
	   //  		g0 = Fold;
	   //  		g1 = Fnew;
	   //  		gprime0 = gradf * dX; // these guys need some work

	   //  		// how to update l1, l2?...


	   //  		l = - gprime0 / (2.0*(g1 - g0 - gprime0));

	   //  		if(l < 0.1)
	   //  			l = 0.1;
	   //  		if(l > 0.5)
	   //  			l = 0.5;
	   //  	}


	   //  }
	   //  else if (f < fold)
	   //  {
	   //  	 l = 1.0;   // use full Newton step
	   //  }
	   //  else if (abs(f-fold) < tolF)
	   //  {
	   //  	if(abs(F) > tolF)
	   //  	{
	   //  		throw "Hit a local min f\n";
	   //  	}
	   //  }












	    // update X_i
	    for(int k = 0; k < n; k++) X[k] += (l * F[k]);

	    // calculate L2-norm of dX iteration
	    dXnorm2 = sqrt(fabs(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]));

		// update individual variables
		lambdai = X[0];
   	 	axi = X[1];
    	azi = X[2];

		i++;

		if(lambdai < 0.0)
		{
			cout << "At iteration " << i << ": ";
    		throw "lambda is negative\n";
		}
    	if(axi < 0.0)
    	{
    		cout << "At iteration " << i << ": ";
    		throw "ax is negative\n";
    	}
    	if(azi < 0.0)
    	{
    		cout << "At iteration " << i << ": ";
    		throw "az is negative\n";
    	}

    	// if(i == Nmax -1)
    	// {
    	// 	cout << ""
    	// 	lambdai = Xold[0];
    	// 	axi = Xold[1];
    	// 	i = 0;
    	// }


	}while((dXnorm2 > toldX) && (Fnorm2 > tolF) && (i < Nmax));



	if(i < Nmax)
	{
		// final answer
		*lambda = lambdai;
		*ax = axi;
		*az = azi;
		cout << i << " ";
	}
	else
	{
		throw "Couldn't find anisotropic variables...\n";
	}

	// free allocated memory

	free_2D(J,n);
}










