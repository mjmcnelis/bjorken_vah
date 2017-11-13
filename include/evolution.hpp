
#include <stdlib.h>

#ifndef EVOLUTION_H

#define EVOLUTION_H

#define ETAS 0.2;  // specific shear viscosity


double dTtt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);

double dTtx_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);

double dTty_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);

double dTtn_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);

double dpl_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);

double dpt_dtau(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double tau);


#endif