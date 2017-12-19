
#include <stdlib.h>
#include <math.h>
#include <string>
#include <string.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include <sstream>
#include <fstream>
#include "evolution.hpp"
#include "qcd.hpp"
#include "anisotropic_transport.hpp"
#include "anisotropic_integrands.hpp"
#include "gauss_integration.hpp"



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//                 ANISOTROPIC TRANSPORT COEFFICIENTS               ::
//                                                                  ::
//     Moments of anisotropic distribution function that appear     ::
//     in transport coefficients of anisotropic hydrodynamics.      ::
//     													            ::
//       				  I240 			I221                        ::
//																	::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


double I240_function(double lambda, double ax, double az, double mbar)
{
		// const int pbar_pts = 16;
		// // a = 2 gauss laguerre roots/weight
		// double pbar_root[pbar_pts] = {0.377613508344741,1.01749195760257,1.94775802042424,3.17692724488987,4.7162400697918,6.58058826577491,8.78946527064707,11.3683230828333,14.350626727437,17.7810957248416,21.7210847965713,26.2581386751111,31.5245960042758,37.7389210025289,45.3185461100898,55.3325835388358};
		// double pbar_weight[pbar_pts] = {0.0486064094670787,0.29334739019044,0.583219363383551,0.581874148596173,0.33818053747379,0.12210596394498,0.0281146258006637,0.00414314919248226,0.000385648533767438,2.20158005631091e-05,7.34236243815652e-07,1.32646044204804e-08,1.15266648290843e-10,3.94706915124609e-13,3.63797825636053e-16,3.45457612313612e-20};

		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

		double pbar_weight[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

		double g = 51.4103536012791;

		double factor_I240 = g * (ax*ax) * (az*az*az*az*az) * (lambda*lambda*lambda*lambda) / (4.0*M_PI*M_PI);

		double answer = factor_I240 * Gauss_Aniso_1D(I240_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		// check conformal formula later:

		return answer;
}



double I221_function(double lambda, double ax, double az, double mbar)
{
		// const int pbar_pts = 16;
		// // a = 2 gauss laguerre roots/weights
		// double pbar_root[pbar_pts] = {0.377613508344741,1.01749195760257,1.94775802042424,3.17692724488987,4.7162400697918,6.58058826577491,8.78946527064707,11.3683230828333,14.350626727437,17.7810957248416,21.7210847965713,26.2581386751111,31.5245960042758,37.7389210025289,45.3185461100898,55.3325835388358};
		// double pbar_weight[pbar_pts] = {0.0486064094670787,0.29334739019044,0.583219363383551,0.581874148596173,0.33818053747379,0.12210596394498,0.0281146258006637,0.00414314919248226,0.000385648533767438,2.20158005631091e-05,7.34236243815652e-07,1.32646044204804e-08,1.15266648290843e-10,3.94706915124609e-13,3.63797825636053e-16,3.45457612313612e-20};

		const int pbar_pts = 32;
		// // a = 2 gauss laguerre roots/weights
		double pbar_root[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

		double pbar_weight[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

		const double g = 51.4103536012791; // degeneracy factor g (nf = 3 flavors)

		double factor_I221 = g * (ax*ax*ax*ax) * (az*az*az) * (lambda*lambda*lambda*lambda) / (8.0*M_PI*M_PI);

		double answer = factor_I221 * Gauss_Aniso_1D(I221_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		// check conformal formula later:

		return answer;
}


double I020_function(double lambda, double ax, double az, double mbar)
{
		const int pbar_pts = 32;

		// a = 0 gauss laguerre roots/weights
		double pbar_root[pbar_pts] = {0.044489365833267,0.234526109519619,0.576884629301886,1.07244875381782,1.72240877644465,2.52833670642579,3.49221327302199,4.61645676974977,5.90395850417424,7.35812673318624,8.9829409242126,10.78301863254,12.7636979867427,14.9311397555226,17.2924543367153,19.8558609403361,22.6308890131968,25.6286360224592,28.8621018163235,32.3466291539647,36.100494805752,40.1457197715394,44.5092079957549,49.2243949873086,54.3337213333969,59.892509162134,65.975377287935,72.6876280906627,80.1874469779135,88.7353404178924,98.829542868284,111.751398097938};

		double pbar_weight[pbar_pts] = {0.109218341952385,0.210443107938813,0.235213229669848,0.195903335972881,0.129983786286072,0.0705786238657174,0.0317609125091751,0.0119182148348386,0.00373881629461152,0.000980803306614955,0.000214864918801364,3.92034196798795e-05,5.93454161286863e-06,7.41640457866755e-07,7.60456787912078e-08,6.35060222662581e-09,4.28138297104093e-10,2.30589949189134e-11,9.79937928872709e-13,3.23780165772927e-14,8.17182344342072e-16,1.54213383339382e-17,2.11979229016362e-19,2.05442967378805e-21,1.3469825866374e-23,5.66129413039736e-26,1.41856054546304e-28,1.91337549445422e-31,1.19224876009822e-34,2.67151121924014e-38,1.33861694210626e-42,4.51053619389897e-48};

		const double g = 51.4103536012791; // degeneracy factor g (nf = 3 flavors)

		double factor_I020 = g * (ax*ax) * (az*az*az) * (lambda*lambda) / (4.0*M_PI*M_PI);

		double answer = factor_I020 * Gauss_Aniso_1D(I020_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		return answer;
}

double I001_function(double lambda, double ax, double az, double mbar)
{
		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.044489365833267,0.234526109519619,0.576884629301886,1.07244875381782,1.72240877644465,2.52833670642579,3.49221327302199,4.61645676974977,5.90395850417424,7.35812673318624,8.9829409242126,10.78301863254,12.7636979867427,14.9311397555226,17.2924543367153,19.8558609403361,22.6308890131968,25.6286360224592,28.8621018163235,32.3466291539647,36.100494805752,40.1457197715394,44.5092079957549,49.2243949873086,54.3337213333969,59.892509162134,65.975377287935,72.6876280906627,80.1874469779135,88.7353404178924,98.829542868284,111.751398097938};

		double pbar_weight[pbar_pts] = {0.109218341952385,0.210443107938813,0.235213229669848,0.195903335972881,0.129983786286072,0.0705786238657174,0.0317609125091751,0.0119182148348386,0.00373881629461152,0.000980803306614955,0.000214864918801364,3.92034196798795e-05,5.93454161286863e-06,7.41640457866755e-07,7.60456787912078e-08,6.35060222662581e-09,4.28138297104093e-10,2.30589949189134e-11,9.79937928872709e-13,3.23780165772927e-14,8.17182344342072e-16,1.54213383339382e-17,2.11979229016362e-19,2.05442967378805e-21,1.3469825866374e-23,5.66129413039736e-26,1.41856054546304e-28,1.91337549445422e-31,1.19224876009822e-34,2.67151121924014e-38,1.33861694210626e-42,4.51053619389897e-48};

		const double g = 51.4103536012791; // degeneracy factor g (nf = 3 flavors)

		double factor_I001 = g * (ax*ax*ax*ax) * (az) * (lambda*lambda) / (8.0*M_PI*M_PI);

		double answer = factor_I001 * Gauss_Aniso_1D(I001_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		return answer;
}


double I000_function(double lambda, double ax, double az, double mbar)
{
		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.044489365833267,0.234526109519619,0.576884629301886,1.07244875381782,1.72240877644465,2.52833670642579,3.49221327302199,4.61645676974977,5.90395850417424,7.35812673318624,8.9829409242126,10.78301863254,12.7636979867427,14.9311397555226,17.2924543367153,19.8558609403361,22.6308890131968,25.6286360224592,28.8621018163235,32.3466291539647,36.100494805752,40.1457197715394,44.5092079957549,49.2243949873086,54.3337213333969,59.892509162134,65.975377287935,72.6876280906627,80.1874469779135,88.7353404178924,98.829542868284,111.751398097938};

		double pbar_weight[pbar_pts] = {0.109218341952385,0.210443107938813,0.235213229669848,0.195903335972881,0.129983786286072,0.0705786238657174,0.0317609125091751,0.0119182148348386,0.00373881629461152,0.000980803306614955,0.000214864918801364,3.92034196798795e-05,5.93454161286863e-06,7.41640457866755e-07,7.60456787912078e-08,6.35060222662581e-09,4.28138297104093e-10,2.30589949189134e-11,9.79937928872709e-13,3.23780165772927e-14,8.17182344342072e-16,1.54213383339382e-17,2.11979229016362e-19,2.05442967378805e-21,1.3469825866374e-23,5.66129413039736e-26,1.41856054546304e-28,1.91337549445422e-31,1.19224876009822e-34,2.67151121924014e-38,1.33861694210626e-42,4.51053619389897e-48};

		const double g = 51.4103536012791; // degeneracy factor g (nf = 3 flavors)

		double factor_I000 = g * (ax*ax) * (az) * (lambda*lambda) / (4.0*M_PI*M_PI);

		double answer = factor_I000 * Gauss_Aniso_1D(I000_integrand, pbar_root, pbar_weight, pbar_pts, ax, az, mbar);

		return answer;
}


double I21_function(double T, double mbar)
{
		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

		double pbar_weight[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

		double g = 51.4103536012791;

		double factor_I21 = g * (T*T*T*T) / (6.0*M_PI*M_PI);

		double answer = factor_I21 * Gauss_Thermal_1D(I21_integrand, pbar_root, pbar_weight, pbar_pts, mbar);

		return answer;
}


double I20_function(double T, double mbar)
{
		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.196943922146667,0.529487866050161,1.01026981913845,1.640616191672,2.42200673335506,3.35625823737525,4.44557319147359,5.69257570606939,7.10035048878373,8.67248915845674,10.413146435518,12.3271087558129,14.4198784243951,16.6977773650005,19.1680758788069,21.839153763432,24.7207039368187,27.823992811746,31.1621978174102,34.7508519173206,38.6084399084037,42.7572156420076,47.2243504952188,52.0435960848824,57.257778984273,62.9227106235616,69.1136582681551,75.9368320953467,83.5517824825995,92.221284870548,102.447989923982,115.52490220024};

		double pbar_weight[pbar_pts] = {0.00825033790777967,0.0671033262747106,0.206386098255352,0.368179392999486,0.446389764546666,0.397211321904435,0.270703020914857,0.144937243765141,0.0619302157291065,0.0213227539141068,0.00594841159169929,0.00134795257769464,0.000248166548996264,3.7053223540482e-05,4.47057760459712e-06,4.33555258401213e-07,3.35571417159735e-08,2.05432200435071e-09,9.83646900727572e-11,3.63364388210833e-12,1.01834576904109e-13,2.12110313498633e-15,3.20100105319804e-17,3.39007439648141e-19,2.41904571899768e-21,1.10270714408855e-23,2.98827103874582e-26,4.34972188455989e-29,2.92108431650778e-32,7.0533942409897e-36,3.81617106981223e-40,1.39864930768275e-45};

		double g = 51.4103536012791;

		double factor_I20 = g * (T*T*T*T) / (2.0*M_PI*M_PI);

		double answer = factor_I20 * Gauss_Thermal_1D(I20_integrand, pbar_root, pbar_weight, pbar_pts, mbar);

		return answer;
}


double I00_function(double T, double mbar)
{
		const int pbar_pts = 32;

		double pbar_root[pbar_pts] = {0.044489365833267,0.234526109519619,0.576884629301886,1.07244875381782,1.72240877644465,2.52833670642579,3.49221327302199,4.61645676974977,5.90395850417424,7.35812673318624,8.9829409242126,10.78301863254,12.7636979867427,14.9311397555226,17.2924543367153,19.8558609403361,22.6308890131968,25.6286360224592,28.8621018163235,32.3466291539647,36.100494805752,40.1457197715394,44.5092079957549,49.2243949873086,54.3337213333969,59.892509162134,65.975377287935,72.6876280906627,80.1874469779135,88.7353404178924,98.829542868284,111.751398097938};

		double pbar_weight[pbar_pts] = {0.109218341952385,0.210443107938813,0.235213229669848,0.195903335972881,0.129983786286072,0.0705786238657174,0.0317609125091751,0.0119182148348386,0.00373881629461152,0.000980803306614955,0.000214864918801364,3.92034196798795e-05,5.93454161286863e-06,7.41640457866755e-07,7.60456787912078e-08,6.35060222662581e-09,4.28138297104093e-10,2.30589949189134e-11,9.79937928872709e-13,3.23780165772927e-14,8.17182344342072e-16,1.54213383339382e-17,2.11979229016362e-19,2.05442967378805e-21,1.3469825866374e-23,5.66129413039736e-26,1.41856054546304e-28,1.91337549445422e-31,1.19224876009822e-34,2.67151121924014e-38,1.33861694210626e-42,4.51053619389897e-48};

		double g = 51.4103536012791;

		double factor_I00 = g * (T*T) / (2.0*M_PI*M_PI);

		double answer = factor_I00 * Gauss_Thermal_1D(I00_integrand, pbar_root, pbar_weight, pbar_pts, mbar);

		return answer;
}





