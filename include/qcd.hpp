
#include <stdlib.h>

#ifndef QCD_H

#define QCD_H

//#define CONFORMAL_EOS
#define BEST

// ideal gas of massless quarks and gluons
#define EOS_FACTOR 15.6269   // Nc=3, Nf=3

//#define CONSTANT_VISCOSITY
#define CONSTANT_ETAS 0.08     // Brookhaven
#define ETAS_MIN 0.08
#define ETAS_SLOPE 0.167728    // Duke
#define ZETA_NORM 1.00
#define T_PEAK (0.155 * 5.067731)

double equilibriumPressure(double e);
double speedOfSoundSquared(double e);
double effectiveTemperature(double e);
double equilibriumEnergyDensity(double T);

double bulkViscosityToEntropyDensity(double T);
double shearViscosityToEntropyDensity(double T);

// quasiparticle thermal functions
double z_Quasiparticle(double T);
double mdmdT_Quasiparticle(double T);
double mdmde_Quasiparticle(double e);
double equilibriumKineticPressure(double T);
double equilibriumBquasi(double T);

// quasiparticle beta function
double beta_shear(double T);
double beta_bulk(double T);


#endif