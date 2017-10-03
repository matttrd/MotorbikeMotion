//
//  roll_angle_formula.c
//  Roll_Angle_Formula [please see the equations in file note.pdf]

#include "roll_angle_formula.h"

#include "roll_angle_formula.h"
#include <math.h>
#include <stdlib.h>
#include "EKFstruct.h"

MOTO *init_moto_default() {
	
    MOTO *ret = (MOTO*) malloc( sizeof(MOTO) );
    
    ret->m = 261.31; //kg
    ret->I_x = 28.088; //kg*m^2
    ret->h1 = 0.6622; // m
    ret->r1 = 0.309; //m
    ret->r2 = 0.287; //m
    ret->r_t = 0.098; //m
    ret->r_r = ret->r1 - ret->r_t;
    ret->Iw1 = 0.5829;// rear, da rivedere
    ret->Iw2 = 0.4065;// front
    ret->Iw = ret->Iw1 + ret->Iw2* pow((ret->r1/ret->r2),2);
    return ret;
}

void create_f(EKF *ekf, double dt) {
	
	//double *out_f = ekf->f;
	
	double phi = ekf->x[0];
	double w_phi = ekf->x[1];
	double u = ekf->u[0];
	double w_z = ekf->u[1];
	
	MOTO *moto = ekf->parameters;
	
	// output f
    ekf->f[0] = w_phi*dt;
    
    double m = moto->m;
    double I_x = moto->I_x;
    double h1 = moto->h1;
    double r_t = moto->r_t;
    double r_r = moto->r_r;
    double I_w = moto->Iw;
    
    
    double inv_inert = 1.0/(I_x + m*h1*(h1 + r_t*cos(phi)) );
    double first_energy = I_w*w_z * (w_z*sin(phi)/cos(phi) - u/(r_r + r_t*(cos(phi) -1)));
    double second_energy = m*(h1*(r_t* pow(w_phi,2)*sin(phi) - w_z*u + 9.81*sin(phi)) - (u*r_t*w_z)/(cos(phi)));
    double third_energy = (m*( pow(w_z,2) )*h1*sin(phi)*(h1*cos(phi) + r_t))/( pow(cos(phi),2) );
    
    ekf->f[1] = inv_inert * (first_energy + second_energy + third_energy)*dt;
	
	// output F
	
	ekf->F[0] = 1;
	ekf->F[1] = dt;
	

	ekf->F[2] = dt*(h1*m*r_t*(I_w*w_z*(-u/(r_r + r_t*(cos(phi) - 1)) + w_z*sin(phi)/cos(phi)) + h1*m*pow(w_z,2)*(h1*cos(phi) + r_t)*sin(phi)/pow(cos(phi),2) + m*(h1*(r_t*pow(w_phi,2)*sin(phi) - u*w_z + 9.81*sin(phi)) - r_t*u*w_z/cos(phi)))*sin(phi)/pow((I_x + h1*m*(h1 + r_t*cos(phi))), 2) + (I_w*w_z*(-r_t*u*sin(phi)/pow((r_r + r_t*(cos(phi) - 1)), 2) + w_z*pow(sin(phi), 2)/pow(cos(phi), 2) + w_z) - pow(h1, 2)*m*pow(w_z, 2)*pow(sin(phi)/cos(phi), 2) + 2*h1*m*pow(w_z, 2)*(h1*cos(phi) + r_t)*pow(sin(phi), 2)/pow(cos(phi), 3) + h1*m*pow(w_z, 2)*(h1*cos(phi) + r_t)/cos(phi) + m*(h1*(r_t*pow(w_phi, 2)*cos(phi) + 9.81*cos(phi)) - r_t*u*w_z*sin(phi)/pow(cos(phi), 2))/(I_x + h1*m*(h1 + r_t*cos(phi)))));
	
	ekf->F[3] = 2*dt*h1*m*r_t*w_phi*sin(phi)/(I_x + h1*m*(h1 + r_t*cos(phi))) + 1.0;
}

void create_h(EKF *ekf) {
    ekf->h[0] = ekf->x[1];

    ekf->H[0] = 0;
    ekf->H[1] = 1;
}

