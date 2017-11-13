
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "anisotropicvariables.hpp"
#include "qcd.hpp"
#include "anisotropic_integrands.hpp"
#include "gauss_integration.hpp"
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
//	   the quantities (e <-> e_kinetic,pl,pt)       		        ::
//                                                                  ::
//                     get_anisotropic_variables                    ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


void get_anisotropic_variables(double e, double pl, double pt, double *lambda, double *ax, double *az)
{
	// this is not designed for conformal mode; it will probably crash
	// could do an ifdef statement



	// (e,pt,pl) should already be updated (not previous values) before running this
	// order: update conserved variables, update inferred variables, update anisotropic variables

	const double T = effectiveTemperature(e);					   // temperature
	const double ekinetic = equilibriumKineticEnergyDensity(T);    // quasiparticle kinetic energy density

	// macropscopic input variables I need to invert
	// are the fa kinetic energy and pressures
	const double Ea = ekinetic,   PTa = pt,	  PLa = pl;

	const double g = 51.4103536012791;                       // degeneracy factor g (nf = 3 flavors)


	const int pbar_pts = 16;
	// a = 2 Gauss Laguerre roots/weights (for F)
	double pbar_rootF[pbar_pts] = {0.377613508344741,1.01749195760257,1.94775802042424,3.17692724488987,4.7162400697918,6.58058826577491,8.78946527064707,11.3683230828333,14.350626727437,17.7810957248416,21.7210847965713,26.2581386751111,31.5245960042758,37.7389210025289,45.3185461100898,55.3325835388358};
	double pbar_weightF[pbar_pts] = {0.0486064094670787,0.29334739019044,0.583219363383551,0.581874148596173,0.33818053747379,0.12210596394498,0.0281146258006637,0.00414314919248226,0.000385648533767438,2.20158005631091e-05,7.34236243815652e-07,1.32646044204804e-08,1.15266648290843e-10,3.94706915124609e-13,3.63797825636053e-16,3.45457612313612e-20};

	// a = 3 gauss laguerre roots/weights (for J)
	double pbar_rootJ[pbar_pts] = {0.567443458991574,1.33290773275989,2.38148248007006,3.72382664209343,5.37212395216187,7.34193662826135,9.65333213726123,12.3323014070182,15.4128500654077,18.9402755758603,22.9765957156019,27.6101814474261,32.9745092032585,39.2898232533641,46.9768962767103,57.1135140237535};
	double pbar_weightJ[pbar_pts] = {0.0650981121009449,0.565273224236402,1.48989313857907,1.86124704489653,1.30047554848272,0.547819646856915,0.143843043208106,0.0237498472421623,0.00244244469350611,0.000152340599558851,5.50126530501356e-06,1.06842629643597e-07,9.92524335897222e-10,3.61863342780077e-12,3.54373278073215e-15,3.58180355287779e-19};




	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




	double lambdai = *lambda;			                     // initial guess for anisotropic variables
	double axi = *ax;										 // are previous time step values
	double azi = *az;


	const int n = 3; 										 // dimension space
	double X[n] = {lambdai, axi, azi};						 // initial guess vector
  	double F[n];  											 // F vector (root equation: F[X] = 0)
	// J = Jacobian of F
	double **J = (double **) malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++) J[i] = (double *) malloc(n* sizeof(double));
 	int pvector[n];									  		 // permutation vector


 	// mbar = m(T) / lambda
	const double thermal_mass = z_Quasiparticle(T) * T;	     // m(T) fixed wpt lambda
	double mbari;                   						 // mbar = m(T) / lambda


	// prefactors and factors of F/J elements: to be evaluated at ith iteration of X
 	double commonfactori;
  	double lambdai2, lambdai3, axi2, azi2, lambdaiaxi3, lambdaiazi3;
  	double factorEai, factorPTai, factorPLai;
  	double factorI2001, factorI2011, factorI2201, factorI401m1,
  	factorI420m1, factorI402m1, factorI421m1, factorI440m1;


	// anisotropic functions evaluated at ith iteration of X
	double Eai, PTai, PLai;
	double I2001, I2011, I2201, I401m1, I420m1, I402m1, I421m1, I440m1;



	int i = 0;				  // starting ith iteration
	int Nmax = 100;			  // max number of iterations
	double dXnorm2;           // L2-norm of dX iteration
	double Fnorm2;		      // L2-norm of F
	double tolX = 1.0e-7;     // tolerance for X
	double tolF = 1.0e-10;    // what's the scale I should use?

	// Find anisotropic variables using 3D Newton Method (change to Broydyn eventually...)
	do{
		// Evaluate mass parameter
		mbari = thermal_mass / lambdai;

		// Evaluate factors and prefactors of F/J
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
	    J[0][2] = I420m1 / lambdaiazi3;
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


	    // update X_i
	    for(int k = 0; k < n; k++) X[k] += F[k];


	    // calculate L2-norm of dX iteration
	    dXnorm2 = sqrt(fabs(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]));


		// update individual variables
		lambdai = X[0];
   	 	axi = X[1];
    	azi = X[2];

		i++;

	} while((dXnorm2 > tolX) && (i < Nmax));


	///cout << i;


	// final answer
	*lambda = lambdai;
	*ax = axi;
	*az = azi;

	// free allocated memory
	free_2D(J,n);
}










