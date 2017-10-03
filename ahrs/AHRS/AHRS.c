//========================================================================================================
//
// Implementation of the quaternion-based indirect-EKF AHRS

========================================================================================================
//
//--------------------------------------------------------------------------------------------------------
// Header files
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// include GSL
//#include </usr/local/include/gsl/gsl_vector.h>
//#include </usr/local/include/gsl/gsl_matrix.h>

//macros
#include "AHRS.h"

//include utility files
#include "AHRSstruct.h"
#include "initialization.h"
#include "mechanization.h"




//--------------------------------------------------------------------------------------------------------
// Variable definitions
//
// Initialization variables
int init_counter = 0;

float nu_g = 0; // Gyro measure variance

// Gravity vector
const float g_n[] = {0.0, 0.0, -1.0};
//Magnetic flux vector
const float m_n[] = {m_n_x, m_n_y, m_n_z};

// AR1 models
//		for gyro
const float sigma_omega_g[] = {sigma_omega_g_x, sigma_omega_g_y, sigma_omega_g_z};
const float F_g[3][3] = {{lambda_g_x, 0.0, 0.0},{0.0, lambda_g_y, 0.0},{0.0, 0.0 ,lambda_g_z}};
const float c_vector[] = {c_vector_x, c_vector_y, c_vector_z};


// Acc validity variables
float mu = INIT_MU; // consider moving it into struct

//---------------------------------------------------------------------------------------------------
// Function declarations

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update

AHRSstatus* const AHRS;


void EKF_AHRS_update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz){
	if (init_counter <= INIT_STOP_COUNTER){
		
		// Initial values for the variables in AHRS
		if (init_counter == 0){
			
			//AHRS = malloc(sizeof(AHRSstatus));
			
			AHRS -> hat_x_g[0] = 0.0;
			AHRS -> hat_x_g[1] = 0.0;
			AHRS -> hat_x_g[2] = 0.0;
			AHRS -> bar_g_b[0] = 0.0;
			AHRS -> bar_g_b[1] = 0.0;
			AHRS -> bar_g_b[2] = 0.0;
			
			AHRS -> bar_m_b[0] = 0.0;
			AHRS -> bar_m_b[1] = 0.0;
			AHRS -> bar_m_b[2] = 0.0;
		}
		
		initialization_stage_estimates(AHRS, gx, gy, gz, ax, ay, az, mx, my, mz, init_counter);
		init_counter++;
		
		if (init_counter == INIT_STOP_COUNTER){
			EKF_initialization(AHRS);
		}
	} else {
		mechanization_step(AHRS, gx, gy, gz, dt);
	}
}




