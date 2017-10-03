//  roll_angle_formula.h

#ifndef roll_angle_formula_h
#define roll_angle_formula_h

#include <stdio.h>
#include "EKFstruct.h"

typedef struct moto_config {
	double m;
	double I_x;
	double h1;
	double r1;
	double r2;
	double r_t;
	double r_r;
	double Iw1;
	double Iw2;
	double Iw;
}MOTO;

MOTO *init_moto_default();

void create_f(EKF *ekf, double dt);

void create_h(EKF *ekf);


#endif /* roll_angle_formula_h */
