
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include "qcd.hpp"


// lattice qcd quantities: expressions are rational polynomial fits to lattice data

double equilibriumPressure(double e) // Peq in units of fm^-4
{
	#ifndef CONFORMAL_EOS
	    // Equation of state from the Wuppertal-Budapest collaboration
	    double e1 = (double)e;
	    double e2 = e*e;
	    double e3 = e2*e;
	    double e4 = e3*e;
	    double e5 = e4*e;
	    double e6 = e5*e;
	    double e7 = e6*e;
	    double e8 = e7*e;
	    double e9 = e8*e;
	    double e10 = e9*e;
	    double e11 = e10*e;
	    double e12 = e11*e;

		double a0 = -0.25181736420168666;
		double a1 = 9737.845799644809;
		double a2 = 1.077580993288114e6;
		double a3 = 3.1729694865420084e6;
		double a4 = 1.6357487344679043e6;
		double a5 = 334334.4309240126;
		double a6 = 41913.439282708554;
		double a7 = 6340.448389300905;
		double a8 = 141.5073484468774;
		double a9 = 0.7158279081255019;
		double a10 = 0.0009417586777847889;
		double a11 = 3.1188455176941583e-7;
		double a12 = 1.9531729608963267e-11;
		double a = (double)fma(a12,e12,fma(a11,e11,fma(a10,e10,fma(a9,e9,fma(a8,e8,fma(a7,e7,fma(a6,e6,fma(a5,e5,fma(a4,e4,fma(a3,e3,fma(a2,e2,fma(a1,e1,a0))))))))))));

		double b0 = 45829.44617893836;
		double b1 = 4.0574329080826794e6;
		double b2 = 2.0931169138134286e7;
		double b3 = 1.3512402226067686e7;
		double b4 = 1.7851642641834426e6;
		double b5 = 278581.2989342773;
		double b6 = 26452.34905933697;
		double b7 = 499.04919730607065;
		double b8 = 2.3405487982094204;
		double b9 = 0.002962497695527404;
		double b10 = 9.601103399348206e-7;
		double b11 = 5.928138360995685e-11;
		double b12 = 3.2581066229887368e-18;
		double b = (double)fma(b12,e12,fma(b11,e11,fma(b10,e10,fma(b9,e9,fma(b8,e8,fma(b7,e7,fma(b6,e6,fma(b5,e5,fma(b4,e4,fma(b3,e3,fma(b2,e2,fma(b1,e1,b0))))))))))));

	    return a / b;
	#else
	    return e / 3.0;
	#endif
}


double speedOfSoundSquared(double e)   // c_s^2 is unitless
{
	#ifndef CONFORMAL_EOS
		// Speed of sound from the Wuppertal-Budapest collaboration
		double e1 = (double)e;
		double e2 = e * e1;
		double e3 = e2 * e1;
		double e4 = e3 * e1;
		double e5 = e4 * e1;
		double e6 = e5 * e1;
		double e7 = e6 * e1;
		double e8 = e7 * e1;
		double e9 = e8 * e1;
		double e10 = e9 * e1;
		double e11 = e10 * e1;
		double e12 = e11 * e1;
		double e13 = e12 * e1;

		return (5.191934309650155e-32 + 4.123605749683891e-23 * e
				+ 3.1955868410879504e-16 * e2 + 1.4170364808063119e-10 * e3
				+ 6.087136671592452e-6 * e4 + 0.02969737949090831 * e5
				+ 15.382615282179595 * e6 + 460.6487249985994 * e7
				+ 1612.4245252438795 * e8 + 275.0492627924299 * e9
				+ 58.60283714484669 * e10 + 6.504847576502024 * e11
				+ 0.03009027913262399 * e12 + 8.189430244031285e-6 * e13)
				/ (1.4637868900982493e-30 + 6.716598285341542e-22 * e
						+ 3.5477700458515908e-15 * e2 + 1.1225580509306008e-9 * e3
						+ 0.00003551782901018317 * e4 + 0.13653226327408863 * e5
						+ 60.85769171450653 * e6 + 1800.5461219450308 * e7
						+ 15190.225535036281 * e8 + 590.2572000057821 * e9
						+ 293.99144775704605 * e10 + 21.461303090563028 * e11
						+ 0.09301685073435291 * e12 + 0.000024810902623582917 * e13);
	#else
		return 1.0 / 3.0;
	#endif
}


double effectiveTemperature(double e)
{
	#ifndef CONFORMAL_EOS
		// Effective temperature from the Wuppertal-Budapest collaboration
		double e1 = (double)e;
		double e2 = e * e1;
		double e3 = e2 * e1;
		double e4 = e3 * e1;
		double e5 = e4 * e1;
		double e6 = e5 * e1;
		double e7 = e6 * e1;
		double e8 = e7 * e1;
		double e9 = e8 * e1;
		double e10 = e9 * e1;
		double e11 = e10 * e1;

		return (1.510073201405604e-29 + 8.014062800678687e-18 * e
				+ 2.4954778310451065e-10 * e2 + 0.000063810382643387 * e3
				+ 0.4873490574161924 * e4 + 207.48582344326206 * e5
				+ 6686.07424325115 * e6 + 14109.766109389702 * e7
				+ 1471.6180520527757 * e8 + 14.055788949565482 * e9
				+ 0.015421252394182246 * e10 + 1.5780479034557783e-6 * e11)
				/ (7.558667139355393e-28 + 1.3686372302041508e-16 * e
						+ 2.998130743142826e-9 * e2 + 0.0005036835870305458 * e3
						+ 2.316902328874072 * e4 + 578.0778724946719 * e5
						+ 11179.193315394154 * e6 + 17965.67607192861 * e7
						+ 1051.0730543534657 * e8 + 5.916312075925817 * e9
						+ 0.003778342768228011 * e10 + 1.8472801679382593e-7 * e11);
	#else
		return powf(e/EOS_FACTOR,0.25);
	#endif
}


double equilibriumEnergyDensity(double T)   // where is this used?
{
	#ifndef CONFORMAL_EOS
		// Effective temperature from the Wuppertal-Budapest collaboration
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		double T13 = T12 * T1;
		double T14 = T13 * T1;
		double T15 = T14 * T1;
		double T16 = T15 * T1;
		double T17 = T16 * T1;
		double T18 = T17 * T1;
		double T19 = T18 * T1;
		double T20 = T19 * T1;
		double T21 = T20 * T1;
		double T22 = T21 * T1;
		double T23 = T22 * T1;

		return (-0.011958188410851651 + 119.89423098138208 * T1
				- 3156.9475699248055 * T2 + 32732.86844374939 * T3
				- 187899.8994764422 * T4 + 712537.3610845465 * T5
				- 1.557049803609345e6 * T6 + 1.4852519861308339e6 * T7
				+ 532132.6079941876 * T8 - 1.963099445042592e6 * T9
				- 4484.44579242679 * T10 + 1.7984228830058286e6 * T11
				+ 119345.25619517374 * T12 - 1.3499773937058165e6 * T13
				- 207838.4995663606 * T14 + 654970.2138652403 * T15
				- 78643.00334616247 * T16 + 40274.00078068926 * T17
				+ 422619.58977657766 * T18 - 409688.07836393174 * T19
				- 62005.75915066359 * T20 + 46788.14270090656 * T21
				+ 40784.330477857235 * T22 - 12589.47744840392 * T23)
				/ (31630.074365558292 - 127100.88940643385 * T1
						+ 173528.1225422275 * T2 - 39403.297956865215 * T3
						- 85582.57873541754 * T4 + 9320.560804233442 * T5
						+ 50882.74198960172 * T6 + 20335.926473421183 * T7
						- 14897.725710713818 * T8 - 23836.484117457 * T9
						- 13726.013896090335 * T10 + 4517.908673107615 * T11
						+ 18056.19917986404 * T12 + 14954.82860467155 * T13
						+ 2569.623976952738 * T14 - 9304.046211514986 * T15
						- 15606.429173842751 * T16 + 8383.710735812094 * T17
						+ 1591.3177623932843 * T18 - 678.748230997762 * T19
						- 33.58687934953277 * T20 + 3.2520554133126285 * T21
						- 0.19647288043440464 * T22 + 0.005443394551264717 * T23);
	#else
		return EOS_FACTOR * powf(T,4.0);
	#endif
}


double derivativeEnergyDensityWithRespectToTemperature(double T)
{
	#ifndef CONFORMAL_EOS
		// Effective temperature from the Wuppertal-Budapest collaboration
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		double T13 = T12 * T1;
		double T14 = T13 * T1;
		double T15 = T14 * T1;
		double T16 = T15 * T1;
		double T17 = T16 * T1;
		double T18 = T17 * T1;
		double T19 = T18 * T1;
		double T20 = T19 * T1;
		double T21 = T20 * T1;
		double T22 = T21 * T1;
		double T23 = T22 * T1;
		return (-0.031494763937498504 + 242.3362205025897*T1 - 7381.284377630215*T2 + 76140.00711480618*T3 - 318626.8893336883*T4 +
	     816969.4526840467*T5 - 1.4077502131456388e6*T6 + 932276.5758216518*T7 + 2.0843919830720688e6*T8 -
	     4.429807180805173e6*T9 - 490343.5182207036*T10 + 8.358430237661866e6*T11 - 5.736749310726519e6*T12 -
	     5.715418080920016e6*T13 + 1.0482485566550283e7*T14 - 5.426507095176461e6*T15 + 185339.33643178307*T16 +
	     1.0547668565460397e6*T17 - 814776.9940396115*T18 + 588800.42598283*T19 - 317439.80902614293*T20 +
	     97410.42727835444*T21 - 15099.361358702017*T22 + 922.756120432997*T23)/
	   (7177.666743649156 - 30194.86618888758*T1 + 35258.44261771079*T2 + 14018.545227919574*T3 - 40048.72302776879*T4 -
	     19144.61212778094*T5 + 43237.317984412286*T6 + 22974.51594985651*T7 - 45584.20063561873*T8 -
	     4546.064537281824*T9 + 37585.605883021555*T10 - 49259.02549065712*T11 + 51255.950072821266*T12 -
	     19777.85808464634*T13 - 17247.826673271855*T14 + 21164.570981481687*T15 - 7534.395375715286*T16 +
	     333.7081011553467*T17 + 337.2157942452745*T18 - 33.57262661808279*T19 - 12.01884142265821*T20 +
	     2.337285690393443*T21 - 0.11495807325727922*T22 + 0.002611670308392503*T23);
	#else
		return 4.0*EOS_FACTOR*powf(T, 3.0);
	#endif
}


double z_Quasiparticle(double T)
{
	// z = m/T (quasiparticle model; nf = 3 flavors)
	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		double T13 = T12 * T1;
		double T14 = T13 * T1;
		double T15 = T14 * T1;
		double T16 = T15 * T1;
		double T17 = T16 * T1;
		double T18 = T17 * T1;
		double T19 = T18 * T1;
		double T20 = T19 * T1;

		return (2.195527549421445e-14 + 5.014273212142939e-11*T1 + 6.8769768080936324e-9*T2 - 1.3090323372462384e-6*T3 + 0.00007503723419734601*T4 -
     	0.0019335423100477788*T5 + 0.0077581370622797526*T6 + 1.0840238709953247*T7 - 40.41237726417249*T8 + 815.5709414412947*T9 -
     	11117.3113417968*T10 + 107720.35990943392*T11 - 727161.0603316505*T12 + 3.12288627176472e6*T13 - 7.240727201149486e6*T14 +
     	8.813651470346985e6*T15 - 3.535759885234317e6*T16 - 5.454281769326236e6*T17 + 1.0415987187648855e7*T18 -
     	8.0189396235499745e6*T19 + 2.6430627991581215e6*T20)/
   	(1.1541703277540495e-16 + 4.122377967105372e-13*T1 + 1.2528690762234088e-10*T2 - 1.1203883703791743e-8*T3 - 3.65428489042621e-8*T4 +
     	0.00004068369176438496*T5 - 0.0024580222716044254*T6 + 0.08433957015198745*T7 - 2.0279246693088906*T8 +
     	36.56366190486074*T9 - 490.05025567116945*T10 + 4477.3308927313055*T11 - 21467.525530376926*T12 - 23929.310531648207*T13 +
     	889228.088912891*T14 - 4.262596231414713e6*T15 + 1.0642287968093066e7*T16 - 1.643826578603126e7*T17 +
     	1.6160269998993594e7*T18 - 9.40765741315934e6*T19 + 2.502796051035953e6*T20);
   	#else
   		return 0.0;
   	#endif
}


double mdmdT_Quasiparticle(double T)
{
	// m * dm/dT (quasiparticle model; nf = 3 flavors)
	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		double T13 = T12 * T1;
		double T14 = T13 * T1;
		double T15 = T14 * T1;
		double T16 = T15 * T1;
		double T17 = T16 * T1;
		double T18 = T17 * T1;
		return (6.816033132910849e-10 + 1.7790915724868891e-6*T1 - 0.00030685573686689566*T2 + 0.023473842070213607*T3 - 1.0732881083537456*T4 +
     	32.832198021861664*T5 - 709.7072077947313*T6 + 11127.167196570503*T7 - 127334.34497884315*T8 + 1.0505151080327253e6*T9 -
     	6.040988782373799e6*T10 + 2.2838438268810038e7*T11 - 5.285418442808765e7*T12 + 7.498857047120465e7*T13 - 6.56141265603832e7*T14 +
     	3.244851452301533e7*T15 - 4.306248519325042e6*T16 - 4.246216963572507e6*T17 + 1.8321297087034914e6*T18)/
   	(1.6642262354933514e-10 + 2.392737679464936e-8*T1 - 5.770833446723635e-6*T2 + 0.00047650968449680815*T3 - 0.02381082180560632*T4 +
     	0.8393191062044852*T5 - 22.059270410143014*T6 + 436.7671455149582*T7 - 6431.104009203983*T8 + 68566.51769617059*T9 -
     	509112.8351052782*T10 + 2.5006513314870317e6*T11 - 7.741180724418022e6*T12 + 1.5432809530797705e7*T13 - 2.021565552418443e7*T14 +
     	1.7198582484739084e7*T15 - 8.803186353145795e6*T16 + 2.1173893883225573e6*T17 - 25629.370279979517*T18);
	#else
		return 0.0;
	#endif
}


// parameterization for specific bulk viscosity
double bulkViscosityToEntropyDensity(double T)
{
	#ifndef CONFORMAL_EOS
		double a0 = -13.45;
		double a1 = 27.55;
		double a2 = -13.77;

		double lambda1 = 0.9;
		double lambda2 = 0.25;
		double lambda3 = 0.9;
		double lambda4 = 0.22;

		double sigma1 = 0.025;
		double sigma2 = 0.13;
		double sigma3 = 0.0025;
		double sigma4 = 0.022;

		double x = T/T_PEAK;

		double result;

		if(x > 1.05)
			result = lambda1*exp(-(x-1.0)/sigma1) + lambda2*exp(-(x-1.0)/sigma2) + 0.001;
		else if(x < 0.995)
			result = lambda3*exp((x-1.0)/sigma3)+ lambda4*exp((x-1.0)/sigma4) + 0.03;
		else
			result = a0 + a1*x + a2*x*x;

		#ifndef CONSTANT_VISCOSITY
			double zeta_norm = ZETA_NORM;
			return zeta_norm * result;
		#else
			return result;
		#endif
	#else
		return 0.0;
	#endif
}


double shearViscosityToEntropyDensity(double T)
{
	#ifndef CONSTANT_VISCOSITY
		double etas_min = ETAS_MIN;
		double etas_slope = ETAS_SLOPE;

		double Tc = 0.154 * 5.067731;

		double ans;

		if(T > Tc)
			ans = etas_min + etas_slope*(T-Tc);
		else if(T <= Tc)
			ans = etas_min;

		return ans;
	#else
		return CONSTANT_ETAS;
	#endif
}



double equilibriumKineticPressure(double T)
{
	// quasiparticle equilbrium kinetic pressure (w/o B(T) background)
	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;

		return (-0.0034061532635185224 + 0.3919654653634235*T1 - 12.47243978413888*T2 + 182.64750463097647*T3 -
     	1488.4344728717017*T4 + 7423.456308051856*T5 - 24215.91600262109*T6 + 53129.172593806084*T7 -
     	77886.57070794867*T8 + 73086.80203597654*T9 - 39349.698083199306*T10 + 8704.509304111943*T11 +
     	524.9236054377599*T12)/
   	(-373.09343601945625 + 2723.149044911457*T1 - 8931.1442216173*T2 + 17720.881793544573*T3 - 23359.75349679735*T4 +
     	20534.437583407194*T5 - 11060.0718544668*T6 + 2789.001986776349*T7 - 7.950258246319416*T8 +
     	19.396345528223684*T9 - 1.8947738765604223*T10 + 0.11466256020804551*T11 - 0.003202799461833753*T12);
   #else
   		return EOS_FACTOR * powf(T,4.0) / 3.0;
	#endif
}

double equilibriumKineticEnergyDensity(double T)
{
	// quasiparticle equilbrium kinetic energy density (w/o B(T) background)
	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;

		return (0.00009491684322244673 - 0.02361379050334422*T1 + 1.3797449684199115*T2 - 33.82886830786682*T3 +
     440.1303217026217*T4 - 3371.5358373636445*T5 + 16005.535555573097*T6 - 49862.71243687168*T7 +
     104848.21092334883*T8 - 148253.7236840333*T9 + 135623.7414817843*T10 - 72796.06845179314*T11 +
     17516.277968272858*T12)/
   (-6.0801761470907785 - 92.02868632549578*T1 + 900.1072873084156*T2 - 3406.550805299837*T3 + 7503.906637596268*T4 -
     10625.464014402616*T5 + 9692.148179782494*T6 - 5257.106451479053*T7 + 1316.2231944400535*T8 -
     12.72184898938659*T9 + 1.3330012931455673*T10 - 0.07856195382532238*T11 + 0.002090630483135667*T12);
   #else
   		return EOS_FACTOR * powf(T,4.0);
	#endif
}

double beta_pi(double T)
{
	// beta_pi = eta / tau_pi
	// temperature dependence of beta_pi = I32/T in quasiparticle kinetic model
	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		return (8.611704298094825e-6 - 0.001405354006113392*T1 + 0.059924402660583354*T2 - 1.1252821226321765*T3 + 11.505088678768702*T4 - 72.86243933500441*T5 +
     	312.67652000568955*T6 - 919.4553300642775*T7 + 1851.880055174282*T8 - 2520.2910195876684*T9 + 2228.4719999229687*T10 -
     	1163.3696941415265*T11 + 275.1649200477449*T12)/
   	(4.814184219646834 - 23.427763291345297*T1 + 78.10907293901035*T2 - 220.9834748507111*T3 + 463.6211365026477*T4 - 658.0962754530466*T5 +
     	604.0491642749903*T6 - 329.15455711633166*T7 + 83.03984882593176*T8 - 0.34262854377950497*T9 + 0.03765047242502595*T10 -
     	0.0024572359354327667*T11 + 0.0000720828838022736*T12);
   	#else
   		return 4.0 * EOS_FACTOR * powf(T,4.0) / 15.0;
   	#endif
}


double beta_Pi(double T)
{
	// beta_Pi = zeta / tau_Pi
	// temperature dependence of beta_Pi = (I32 - I31*I31/I30)/T in quasiparticle kinetic model
	// the moments are calculated with a kinetic distribution with mass m(T), no fit to lattice

	// I need to review this calculation


	#ifndef CONFORMAL_EOS
		double T1 = (double)T;
		double T2 = T1 * T1;
		double T3 = T2 * T1;
		double T4 = T3 * T1;
		double T5 = T4 * T1;
		double T6 = T5 * T1;
		double T7 = T6 * T1;
		double T8 = T7 * T1;
		double T9 = T8 * T1;
		double T10 = T9 * T1;
		double T11 = T10 * T1;
		double T12 = T11 * T1;
		double T13 = T12 * T1;
		double T14 = T13 * T1;
		// 5/3*betapi - cs2kinetic(e+p)
		return (0.000010905760709167861 - 0.004683989457592804*T1 + 0.24754406178287106*T2 - 4.885445624666566*T3 + 50.09614745592617*T4 -
	     307.7032088120906*T5 + 1199.345936356805*T6 - 3134.7810945740785*T7 + 5705.875248385415*T8 -
	     7386.615125073046*T9 + 6838.523042759372*T10 - 4474.354253326277*T11 + 1997.5827901556845*T12 -
	     562.5161507019361*T13 + 79.67974066311064*T14)/
	   (-40.624628025755804 + 107.6723401385496*T1 - 23.19390975206659*T2 - 95.34803927008286*T3 - 47.579669340536825*T4 +
	     66.02920940835048*T5 + 194.0662782679704*T6 + 14.784256174958234*T7 - 637.8693894211432*T8 +
	     806.4328841485209*T9 - 469.0986465043204*T10 + 137.20558809601678*T11 - 9.98223482226815*T12 +
	     0.30712590209276874*T13 - 0.002441723962961617*T14);

		// update: try the formula expressed in paper
		// return (0.0002245190537071987 + 0.001589813465310276*T1 - 0.032887293474307244*T2 + 0.01218795644310531*T3 + 0.7910874363202375*T4 -
  //    3.6392267208984035*T5 + 7.477646242067765*T6 - 8.532303699790662*T7 + 5.785883636049677*T8 - 2.3663215208386874*T9 +
  //    0.5724763080969358*T10 - 0.0745193511774874*T11 + 0.00407583650278803*T12)/
  //  (-50.32078956062687 + 326.03646162349486*T1 - 934.0564150931393*T2 + 1558.029596953221*T3 - 1680.0755778146636*T4 +
  //    1231.5163043799407*T5 - 628.8439282106788*T6 + 225.5054194952194*T7 - 56.414571553536*T8 + 9.618849291778407*T9 -
  //    1.0621740126341475*T10 + 0.06832115730088792*T11 - 0.0019387124710670907*T12);


		//without mass dependent term: beta_Pi = (I32 - I31*I31/I30)/T
		// return (-2.119949615708305e-6 - 0.000026753813875308168*T + 0.009837838109661427*T2 - 0.27333268346935624*T3 + 3.0296961640990436*T4 -
  //    	16.873010802047027*T5 + 57.05168278130282*T6 - 128.05466864079622*T7 + 201.44986207473065*T8 - 228.2823602875621*T9 +
  //    	183.40857129893806*T10 - 94.54476713724578*T11 + 23.393699268595284*T12)/
  //  	(27.724812013804105 - 265.08864856909366*T + 1176.4992316219568*T2 - 3080.0694129542835*T3 + 5158.11613454253*T4 - 5651.242750961873*T5 +
  //    	3951.553670894849*T6 - 1587.3097532527718*T7 + 243.67506179297166*T8 + 30.902785322903576*T9 - 3.1510577959848893*T10 +
  //    	0.19411208556312864*T11 - 0.005443917177925473*T12);

   		// includes mass dependent term: beta_Pi = (I32 - I31*I31/I30)/T + mdmdT*I11*I31/I30
   	// 	return (-0.001436141900836525 + 0.1260092369163144*T1 - 2.765730686955577*T2 + 26.61638920834812*T3 - 172.63339179115172*T4 + 923.8062244678376*T5 -
    //  	2645.3073324819816*T6 + 2910.434246047423*T7 + 2748.89708956389*T8 - 12032.409752869014*T9 + 14860.25153744898*T10 -
    //  	8647.237394047126*T11 + 2032.1014938266533*T12)/
   	// (1080.0320775870898 - 6099.26859611491*T1               + 15962.390547209012*T2 - 27636.169853785006*T3 + 37894.33912748253*T4 - 41629.88872423066*T5 +
    //  	32426.735100396607*T6 - 14831.011024155318*T7 + 2485.461652649497*T8 + 417.7637232720538*T9 - 47.02302553146102*T10 +
    //  	2.8155670946603135*T11 - 0.07742531859742523*T12);
   	#else
   		return 0.0;
   	#endif
}












