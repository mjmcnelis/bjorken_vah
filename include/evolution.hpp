
#include <stdlib.h>

#ifndef EVOLUTION_H

#define EVOLUTION_H


double dTtt_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double dTtx_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double dTty_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double dTtn_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double dpl_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double dpt_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);

double db_dt(double Ttt, double Ttx, double Tty, double Ttn, double pl, double pt, double b, double ut, double ux, double uy, double un, double e, double p, double lambda, double ax, double az, double t);


#endif