
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
	// takes A and decomposes it into LU (with row permutations P)
	// A = n x n matrix; function does A -> PA = LU; (L,U) of PA stored in same ** array
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



void calcF(double Ea, double PTa, double PLa, double thermal_mass, double X[], double F[], int n)
{
	if(n != 3) throw "check variable dimension\n";

	double lambda = X[0];
	double ax = X[1];
	double az = X[2];

	// guass laguerre roots and weights
	const int pbar_pts = 64;

	// double pbar_rootF[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

	// double pbar_weightF[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};


	double pbar_rootF[pbar_pts] = {0.100682305413897,0.270520836483334,0.51569963592892,0.836488008915462,1.23309742130066,1.70576669355907,2.25477551802657,2.88044846967868,3.58315705712956,4.36332133108097,5.22141146098448,6.15794942036906,7.17351084198293,8.26872707803859,9.44428749245818,10.7009420101719,12.0395039495236,13.4608531663356,14.9659395416811,16.5557868497556,18.2314970474232,19.9942550331282,21.8453339300473,23.7861009568067,25.8180239590467,27.9426786869022,30.1617569174568,32.4770755379225,34.8905867252913,37.4043893822821,40.0207420185323,42.7420773014151,45.5710185441942,48.5103984525402,51.5632805164228,54.7329835166212,58.0231097182629,61.4375774542106,64.9806589672541,68.657024593526,72.4717946462972,76.4306007215915,80.5396586263206,84.8058557708515,89.2368567364209,93.8412319199621,98.6286158193952,103.609903872579,108.797500148616,114.205633163617,119.850764560272,125.752126880183,131.932444832859,138.418924137688,145.244642251931,152.450564030877,160.088570135959,168.226210922997,176.954588288189,186.40237480281,196.76320677131,208.356839722359,221.797108929038,238.670543959774};

	double pbar_weightF[pbar_pts] = {0.00121325475254743,0.0115830307874021,0.0449300036230642,0.108723355347864,0.192554390775891,0.269945158837149,0.313196888807452,0.3091769301254,0.264494623617809,0.198598314986153,0.132082681108754,0.0783314174790193,0.0416307828252456,0.0199027620462326,0.00858345076161478,0.00334641930025591,0.00118125650134425,0.000377947187293001,0.000109686881563897,2.88861374224277e-05,6.90381926262571e-06,1.49727487397793e-06,2.94561519553122e-07,5.25379843675339e-08,8.48925787486123e-09,1.24153688239005e-09,1.64156509022329e-10,1.95975096072436e-11,2.10932526722343e-12,2.04342941412443e-13,1.77843322251126e-14,1.38763902259188e-15,9.6845877796311e-17,6.03052529286172e-18,3.34112985700812e-19,1.64199732576702e-20,7.13411087480561e-22,2.73023949071211e-23,9.16648011305265e-25,2.68792443242171e-26,6.85043816923282e-28,1.50923804767702e-29,2.85714918258405e-31,4.61690767984817e-33,6.32098857645367e-35,7.27139225952023e-37,6.96267213431972e-39,5.49095306851804e-41,3.52351639625593e-43,1.81437534698469e-45,7.37759880699764e-48,2.32475267398332e-50,5.55221674798189e-53,9.78674399903822e-56,1.23270861980389e-58,1.06602396324181e-61,6.01707087138138e-65,2.07559338517285e-68,4.004692345493e-72,3.8146451076973e-76,1.48753138733599e-80,1.74539091081978e-85,3.42627051390618e-91,2.51856644360717e-98};


	double g = 51.4103536012791;            // degeneracy factor
	double mbar = thermal_mass / lambda;    // mbar = m(T) / lambda
  	double ax2 = ax * ax;
  	double az2 = az * az;
  	double commonfactor = g * ax2 * az * lambda * lambda * lambda * lambda / (4.0*M_PI*M_PI);

	double factorEa = commonfactor;
    double factorPTa = 0.5 * commonfactor * ax2;
    double factorPLa = commonfactor * az2;

    double Eai = factorEa * Gauss_Aniso_1D(Ea_integrand, pbar_rootF, pbar_weightF, pbar_pts, ax, az, mbar);
	double PTai = factorPTa * Gauss_Aniso_1D(PTa_integrand, pbar_rootF, pbar_weightF, pbar_pts, ax, az, mbar);
	double PLai = factorPLa * Gauss_Aniso_1D(PLa_integrand, pbar_rootF, pbar_weightF, pbar_pts, ax, az, mbar);

	F[0] = Eai - Ea;
	F[1] = PTai - PTa;
	F[2] = PLai - PLa;
}

void calcJ(double Ea, double PTa, double PLa, double thermal_mass, double X[], double dX[], double F[], double Fcurrent[], double ** J, double ** Jcurrent, jacobian jtype, int n)
{
	if(n != 3) throw "check variable dimension\n";

    switch(jtype)
    {
    	case newton:
    	{
    		// Newton Jacobian
			double lambda = X[0];
			double ax = X[1];
			double az = X[2];

			const int pbar_pts = 64;

			// double pbar_rootJ[pbar_pts] = {0.299618729049241,0.701981065353977,1.24974569814569,1.94514382443706,2.78994869155499,3.78614305529879,4.93605293430173,6.24240692441721,7.70838362900271,9.3376617765878,11.1344784990019,13.1036993239838,15.2509033992287,17.5824881879873,20.1057991514724,22.8292918399089,25.7627366009885,28.9174802233293,32.3067850039226,35.9462752181738,39.8545359681502,44.0539338428991,48.5717701879893,53.4419507581545,58.7074908793654,64.4244418290391,70.6683893377061,77.5459900633602,85.2174664086695,93.9467116599065,104.238552969691,117.391742318923};

			// double pbar_weightJ[pbar_pts] = {0.00660146448073508,0.0813584931042281,0.347537436309438,0.809963198105261,1.22739584119905,1.32050782861975,1.06049919505728,0.655616488144915,0.318173017008472,0.122743109012855,0.0379333897858022,0.00943187028689987,0.00188978713293874,0.000304914974586437,3.95130877631855e-05,4.09377958251348e-06,3.36921618654073e-07,2.1841295448875e-08,1.10337736506627e-09,4.28638379146177e-11,1.25966453444067e-12,2.74423030367617e-14,4.32175197361363e-16,4.76686817705967e-18,3.53643350342934e-20,1.67355018349782e-22,4.70254099995936e-25,7.09116556196869e-28,4.93082516196282e-31,1.23284946609868e-34,6.91389702736573e-39,2.63586492716958e-44};

			double pbar_rootJ[pbar_pts] = {0.154232997940043,0.361092175674186,0.642219428767814,0.998308006794363,1.42968128832083,1.936627238075,2.51945540506056,3.17851259873621,3.91418893108869,4.72692103736593,5.61719439252839,6.58554537429053,7.63256333057594,8.75889277039389,9.96523574365244,11.2523544540111,12.621074140785,14.0722862639348,15.6069520273647,17.2261062788428,18.9308618293541,20.7224142404084,22.6020471347841,24.5711380944953,26.6311652196597,28.7837144336854,31.0304876341733,33.3733118056311,35.8141492301113,38.3551089559958,40.9984597143302,43.7466445076039,46.6022971392952,49.5682610059129,52.6476105394023,55.8436757701773,59.1600705844319,62.6007253800711,66.1699249920999,69.8723529722297,73.713143584806,77.6979432442361,81.8329835994095,86.1251691132362,90.582182855848,95.2126154248177,100.026123569947,105.03362745559,110.2475588866,115.682177808651,121.353981877887,127.282245410662,133.489742240704,140.003736751727,146.857377706354,154.0917184313,161.758752072877,169.926176314968,178.685293322514,188.165061333651,198.559551698188,210.189248551775,223.669409375799,240.590314300395};

			double pbar_weightJ[pbar_pts] = {0.000535770785642018,0.00799959651061838,0.0443923689573271,0.144338079765975,0.328163505853751,0.570603859220408,0.799361263929005,0.932566304404982,0.926461049686763,0.796162953689164,0.59863161158653,0.397181950261682,0.234037500925501,0.123079914167528,0.0579885822872833,0.0245482004819927,0.00935808941944307,0.00321789158123392,0.000999306700534512,0.000280496517163487,7.11984810253681e-05,1.63461238029384e-05,3.39415816843481e-06,6.37222551348658e-07,1.08109474528694e-07,1.65626551320671e-08,2.28918634762506e-09,2.8511605746787e-10,3.19572238503125e-11,3.21849467339972e-12,2.90744971155867e-13,2.35121753664469e-14,1.69840636333814e-15,1.09320366120742e-16,6.25319153711965e-18,3.16922540070663e-19,1.41851235223462e-20,5.58695386460294e-22,1.92864703562425e-23,5.80981448943726e-25,1.5198503239326e-26,3.43432400582798e-28,6.66348873291895e-30,1.10283030827861e-31,1.54544384614492e-33,1.81859582540058e-35,1.78033499283807e-37,1.43467694211503e-39,9.40271157290723e-42,4.94288514535937e-44,2.0510089868996e-46,6.59276910303509e-49,1.60565836109471e-51,2.88533628831282e-54,3.70413330369449e-57,3.26423147162308e-60,1.87730204229254e-63,6.59792136630667e-67,1.29710111095039e-70,1.25916246862726e-74,5.00596106810534e-79,5.99270431933009e-84,1.20186591406375e-89,9.05378208966149e-97};

			double g = 51.4103536012791;
			double mbar = thermal_mass / lambda;

		  	double lambda2 =  lambda * lambda;
		  	double lambda3 = lambda2 * lambda;
		  	double ax2 = ax * ax;
		  	double az2 = az * az;
		  	double lambdaax3 = lambda * ax2 * ax;
		  	double lambdaaz3 = lambda * az2 * az;

		  	double commonfactor = g * ax2 * az * lambda*lambda*lambda*lambda*lambda / (4.0*M_PI*M_PI);

    		// calculate Eai,PTai,PLai from F
			double Eai = F[0] + Ea;
			double PTai = F[1] + PTa;
			double PLai = F[2] + PLa;

	    	double factorI2001 = commonfactor;
		    double factorI2011 = 0.5 * commonfactor * ax2;
			double factorI2201 = commonfactor * az2;
			double factorI402m1 = 0.125 * commonfactor * ax2 * ax2;
			double factorI421m1 = 0.5 * commonfactor * ax2 * az2;
			double factorI440m1 = commonfactor * az2 * az2;

			double I2001 = factorI2001 * Gauss_Aniso_1D(I2001_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    double I2011 = factorI2011 * Gauss_Aniso_1D(I2011_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    double I2201 = factorI2201 * Gauss_Aniso_1D(I2201_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    double I402m1 = factorI402m1 * Gauss_Aniso_1D(I402m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    double I421m1 = factorI421m1 * Gauss_Aniso_1D(I421m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    double I440m1 = factorI440m1 * Gauss_Aniso_1D(I440m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);

		    J[0][0] = I2001/lambda2;	  J[1][0] = I2011/lambda2;		   J[2][0] = I2201/lambda2;
		    J[0][1] = 2.0*(Eai+PTai)/ax;  J[1][1] = 4.0*I402m1/lambdaax3;  J[2][1] = 2.0*I421m1/lambdaax3;
		    J[0][2] = (Eai+PLai)/az;	  J[1][2] = I421m1/lambdaaz3;      J[2][2] = I440m1 / lambdaaz3;

		    // need to change common factor
		    // double factorI401m1 = factorI2011;
			// double factorI420m1 = factorI2201;
			// double I401m1 = factorI401m1 * Gauss_Aniso_1D(I401m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    // double I420m1 = factorI420m1 * Gauss_Aniso_1D(I420m1_integrand, pbar_rootJ, pbar_weightJ, pbar_pts, ax, az, mbar);
		    //J[0][2] = I420m1 / lambdaaz3;
		    //J[0][1] = 2.0 * I401m1 / lambdaiaxi3;

    		break;
    	}
    	case broyden:
    	{
    		// Broyden Jacobian
    		double dXabs = sqrt(dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2]);
    		double dX_unit[3] = {0.0,0.0,0.0};
    		double Jcurrent_dX[3] = {0.0,0.0,0.0};

    		for(int i = 0; i < n; i++)
			{
				dX_unit[i] = dX[i] / dXabs;
				for(int k = 0; k < n; k++) Jcurrent_dX[i] += Jcurrent[i][k]*dX[k];
			}
    		for(int i = 0; i < n; i++)
    			for(int j = 0; j < n; j++)
    				J[i][j] = Jcurrent[i][j] + (F[i] - Fcurrent[i] - Jcurrent_dX[i]) * dX_unit[j] / dXabs;
    				// do I get any numerical error from this?
    		break;
    	}
    	default:
    		throw "specify method for J\n";
    }
}



void linebacktrack(double *l, double Ea, double PTa, double PLa, double thermal_mass, double Xcurrent[], double dX[], double g0, double F[], double gradf[], double toldX, int n)
{
	double ls = 1.0;         // default value
	double alpha = 0.0001;   // descent parameter

	double dXabs = sqrt(dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2]);

	// update F at least once
	double X[3];
	for(int k = 0; k < n; k++)
	{
		X[k] = Xcurrent[k] + ls*dX[k];
	}
	calcF(Ea, PTa, PLa, thermal_mass, X, F, n);


	if(ls < (toldX/dXabs))
	{
		*l = ls;
		return;
	}

	double f = 0.5*(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
	double g1 = f;

	double gprime0 = 0.0;
	for(int k = 0; k < n; k++)
	{
		gprime0 += gradf[k]*dX[k];
	}
	if(gprime0 >= 0.0)
	{
		throw "\nRoundoff issue with gradient descent";
	}

	double lroot, lprev, fprev;
	double a, b, z;
	double N = 10;

	// Main loop
	for(int i = 0; i < N; i++)
	{
		// check for sufficient decrease
		if(f <= g0 + ls*alpha*gprime0)
		{
			*l = ls;
			return;
		}
		// quadratic formula
		if(i == 0)
		{
			lroot = - gprime0 / (2.0*(g1 - g0 - gprime0));
		}
		// cubic formula
		else
		{
			a = ((g1-gprime0*ls-g0)/(ls*ls) - (fprev-gprime0*lprev-g0)/(lprev*lprev)) / (ls-lprev);
			b = (-lprev*(g1-gprime0*ls-g0)/(ls*ls) + ls*(fprev-gprime0*lprev-g0)/(lprev*lprev)) / (ls-lprev);

			if(a == 0.0) // Solve dg/dl = 0 (a = 0)
				lroot = - gprime0 / (2.0*b);
			else
			{
				z = b*b - 3.0*a*gprime0;

				if(z < 0.0)
					lroot = 0.5*ls;
				else
					lroot = (-b + sqrt(z)) / (3.0*a);

				// does result change if I used general formula?

				// if(z < 0.0)
				// 	lroot = 0.5*ls;
				// else if(b <= 0.0)
				// 	lroot = (-b + sqrt(z)) / (3.0*a);
				// else
				// 	lroot = - gprime0 / (b + sqrt(z));  // ?????
			}
			if(lroot > 0.5*ls) lroot = 0.5*ls;

		}


		// store current values for next iteration
		lprev = ls;
		fprev = f;


		// update l and f
		ls = fmax(lroot, 0.1*ls);

		for(int k = 0; k < n; k++)
		{
			X[k] = Xcurrent[k] + ls*dX[k];
		}
		calcF(Ea, PTa, PLa, thermal_mass, X, F, n);
		f = 0.5*(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
	}

	*l = ls;
	return;
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                       ANISOTROPIC VARIAbLES                      ::
//                                                                  ::
//     Compute the anisotropic variables (lambda,ax,az) from        ::
//	   the quantities (e,b,pl,pt)       		                    ::
//                                                                  ::
//                     get_anisotropic_variables                    ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::






void get_anisotropic_variables(double e, double pl, double pt, double b, double *lambda, double *ax, double *az)
{
	// not designed for conformal mode, it will crash b/c detJ = 0
	// (e,b,pt,pl) should already be updated before running this
	// order: update conserved variables, update inferred variables, update anisotropic variables

	const double T = effectiveTemperature(e);				 // temperature
	double thermal_mass = z_Quasiparticle(T) * T;	   	     // quasiparticle mass

	const double Ea = e - b;								 // kinetic energy density
	const double PTa = pt + b;							     // kinetic transverse pressure
	const double PLa = pl + b;							     // kinetic longitudinal pressure

	if(Ea < 0.0) throw "Ea is out of bounds!\n";
	if(PTa < 0.0) throw "PTa is out of bounds!\n";
	if(PLa < 0.0) throw "PLa is out of bounds!\n";

	const int n = 3; 										 // multivariable dimension
	double X[n] = {*lambda, *ax, *az};						 // initialize solution vector to guess; will iterate w/ + l*dX
	double Xcurrent[n] = {*lambda, *ax, *az};			     // holder for current X solution
	double dX[n];							 				 // dX iteration
  	double F[n];  											 // F(X) vector in root equation F(Xsol) = 0
  	double Fcurrent[n];  									 // F(Xcurrent)
  	double f;		                                  	     // f = 0.5 F * F at X
  	double fcurrent;										 // f = 0.5 F * F at Xcurrent
	double **J = (double**)malloc(n*sizeof(double*));        // J = Jacobian of F at X
	double **Jcurrent = (double**)malloc(n*sizeof(double*)); // Jacobian at Xcurrent
	double gradf[n];										 // gradient of f = F*F/2 at Xcurrent
	double gradf_sum; 										 // summation variable for gradf calculation
 	int pvector[n];									  		 // permutation vector for LUP solver
 	jacobian jtype = newton;                                 // jacobian method variable (newton for i = 0)

 	// allocate Jacobian matrix
 	for(int k = 0; k < n; k++)
	{
		J[k] = (double*)malloc(n*sizeof(double));
 		Jcurrent[k] = (double*)malloc(n*sizeof(double));
 	}

	int Nmax = 2000;	      // max number of iterations
	double dXabs;             // magnitude of dX
	double Fabs;		      // magnitude of F
	double toldX = 1.0e-6;    // tolerance for dX ~ 0
	double tolF = 1.0e-12;    // tolerance for F ~ 0
	double tolmin = 1.0e-6;   // tolerance for spurious convergence to local min of f = F*F/2
	double stepmax = 100.0;   // scaled maximum step length allowed in line searches (I don't know what this means...)
	double gprime0;           // descent derivative at Xcurrent
	double l = 1.0; 		  // partial step parameter


	//stepmax calculation (no idea what this means, so far it doesn't seem to be used..)
	stepmax = stepmax*fmax(sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]),sqrt(3.0)*(double)n);


	// compute F at Xcurrent
	calcF(Ea, PTa, PLa, thermal_mass, Xcurrent, F, n);


	// 3D Newton method with line backtracking (main loop)
	for(int i = 0; i < Nmax; i++)
	{
		// run at least once before
		// checking for convergence
		if(i > 0)
		{
			if(dXabs <= toldX & Fabs <= tolF)
			{
				//cout << i << " ";

				// update anisotropic variables
				*lambda = X[0];
				*ax = X[1];
				*az = X[2];

				// free allocated memory
				free_2D(J,n);
				free_2D(Jcurrent,n);

				return;
			}
			else
			{
				jtype = broyden;
			}
		}


	    // compute default J and f at Xcurrent
	    calcJ(Ea, PTa, PLa, thermal_mass, Xcurrent, dX, F, Fcurrent, J, Jcurrent, jtype, n);
	    fcurrent = 0.5*(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);


    	for(int k = 0; k < n; k++)
    	{
    		// store current F and default J
    		Fcurrent[k] = F[k];
    		for(int j = 0; j < n; j++)
    		{
    			Jcurrent[k][j] = J[k][j];
    		}
    		// gradient of f at Xcurrent (gradf_k = F_m * J_mk)
    		gradf_sum = 0.0; // clear sum
	    	for(int m = 0; m < n; m++) gradf_sum += F[m] * J[m][k];
	    	gradf[k] = gradf_sum;
    	}


	    // Solve matrix equation: J * dX = - F
	    for(int k = 0; k < n; k++) F[k] = - F[k];  // change sign of F first
	    // LU solver routine
	    LUP_decomposition(J, n, pvector);          // LUP of J now stored in J (pvector also updated)
	    LUP_solve(J, n, pvector, F);               // F now stores dX
	    for(int k = 0; k < n; k++) dX[k] = F[k];   // load full Newton/Broyden iteration dX


	    // rescale dX if too large
	    dXabs = sqrt(dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2]);
		if(dXabs > stepmax)
		{
			for(int k = 0; k < n; k++) dX[k] *= (stepmax/dXabs);
		}


		// descent derivative
		gprime0 = 0.0;
		for(int k = 0; k < n; k++)
		{
			gprime0 += gradf[k]*dX[k];
		}
		switch(jtype)
		{
			case newton:
			{
				if(gprime0 >= 0.0) throw "error in newton gradient descent\n";
			}
			case broyden:
			{
				break;
			}
			default:
			{
				throw "correct gradient descent\n";
			}
		}



		switch(jtype)
		{
			case newton:
			{
				// line backtracking algorithm (solve for l and update F(Xcurrent+l*dX))
				linebacktrack(&l, Ea, PTa, PLa, thermal_mass, Xcurrent, dX, fcurrent, F, gradf, toldX, n);

				// update X
			    for(int k = 0; k < n; k++)
			    {
			    	X[k] = Xcurrent[k] + l*dX[k];
			    }
				break;
			}
			case broyden:
			{
				// default Broyden iteration (l = 1)
				for(int k = 0; k < n; k++)
				{
					X[k] = Xcurrent[k] + dX[k];
				}
				// F vector and magnitude at default X
				calcF(Ea, PTa, PLa, thermal_mass, X, F, n);
				f = 0.5*(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);


				// If default Broyden step doesn't decrease f,
				// redo the default Newton step with exact Jacobian
				if(f > fcurrent - 0.2*fabs(gprime0))
				{
					// reset F
			    	for(int k = 0; k < n; k++)
			    	{
			    		F[k] = Fcurrent[k];
			    	}

			    	jtype = newton;

			    	// reinitialize J as exact
			    	calcJ(Ea, PTa, PLa, thermal_mass, Xcurrent, dX, F, Fcurrent, J, Jcurrent, jtype, n);


			    	// store J and calculate gradf
			    	for(int k = 0; k < n; k++)
			    	{
			    		for(int j = 0; j < n; j++)
			    		{
			    			Jcurrent[k][j] = J[k][j];
			    		}
				    	gradf_sum = 0.0; // clear sum
				    	for(int m = 0; m < n; m++) gradf_sum += F[m] * J[m][k];
				    	gradf[k] = gradf_sum;
				    }


				    // resolve matrix equation
				    for(int k = 0; k < n; k++) F[k] = - F[k];
				    LUP_decomposition(J, n, pvector);
				    LUP_solve(J, n, pvector, F);
				    for(int k = 0; k < n; k++) dX[k] = F[k];


				    dXabs = sqrt(dX[0]*dX[0]+dX[1]*dX[1]+dX[2]*dX[2]);
					if(dXabs > stepmax)
					{
						for(int k = 0; k < n; k++) dX[k] *= (stepmax/dXabs);
					}


					// default Newton iteration
					for(int k = 0; k < n; k++) X[k] = Xcurrent[k] + dX[k];


					// line backtracking
					linebacktrack(&l, Ea, PTa, PLa, thermal_mass, Xcurrent, dX, fcurrent, F, gradf, toldX, n);


					// redo update for X
				    for(int k = 0; k < n; k++)
				    {
					    X[k] = Xcurrent[k] + l*dX[k];
				    }
				}
				break;
			}
			default:
				throw "correct line backtracking routine\n";
		}


	    // update convergence values
	    Fabs = sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
	    dXabs = sqrt(dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2]);


	    // store current solution
		for(int k = 0; k < n; k++)
		{
			Xcurrent[k] = X[k];
		}


	    // for some reason the solution can yield
	    // negative values at initialization

		if(X[0] < 0.0)
		{
			//cout << "At iteration " << i << ": ";
    		throw "lambda is negative\n";
		}
    	if(X[1] < 0.0)
    	{
    		//cout << "At iteration " << i << ": ";
    		throw "ax is negative\n";
    	}
    	if(X[2] < 0.0)
    	{
    		//cout << "At iteration " << i << ": ";
    		throw "az is negative\n";
    	}

    	//cout << X[0] << endl;

    	// to do: what should I do if I hit a local minimum? (make a |dF| convergence criteria)
    	if(i == Nmax - 1)
		{
			// free allocated memory and throw error
			free_2D(J,n);
			free_2D(Jcurrent,n);
			throw "Iterations exceeded: couldn't find anisotropic variables...\n";
		}
    }
}










