//
//  driver_functions.c
//  GNSS_INS
//
#include <Accelerate/Accelerate.h>

#include "driver_functions.h"
#include "/init/initialization.h"

// global variables for initialization

int init_counter = 0;

#define SPEED_IDX 3
#define COURSE_IDX 4


void updateGNSS_INS(precision *measurements, EKF *ekf, precision dt,GNSS_INS_t* ins_struct){
	
	if (init_counter <= INIT_STOP_COUNTER){
		//precision x0[] = {45.4170, 11.8669, 172.0583, 0, 0, 0, 1, 0, 0, 0, 0.0038, 0.0061, -0.0219};
		
		// Initial values for the variables in the struct EKF

		initialization_stage_estimates(ekf, measurements, init_counter);
		init_counter++;
		
		if (init_counter == INIT_STOP_COUNTER){
			GNSS_INS_initialization(ekf,ins_struct);
		}
	} else {
		// ekf prediction
		ekf_prediction(ekf, dt);
		normalize_quaternion(ekf);
	// copy the first 6 measures
#if PRECISION == DOUBLE
	cblas_dcopy(6, measurements, 1, ekf->z, 1);
#else
	cblas_scopy(6, measurements, 1, ekf->z, 1);
#endif
	// reset measurements in navigation frame
	ekf->z[8] = - ge;

	
	// copy new inputs in ekf->u
#if PRECISION == DOUBLE
	cblas_dcopy(N_INPUTS, measurements, 1, ekf->u, 1);
#else
	cblas_scopy(N_INPUTS, measurements, 1, ekf->u, 1);
#endif
		// ekf update
		ekf_update(ekf);
		normalize_quaternion(ekf);
	}
}


void preprocess_GPS_data(double *GPS_raw_data,GNSS_INS_t *ins_master_struc){
    
    
    double speed  = GPS_raw_data[SPEED_IDX]*KNOTS_TO_METERS;
    double course = GPS_raw_data[COURSE_IDX];
    double v_n = speed*cos(course*M_PI/180.0);
    double v_e = speed*sin(course*M_PI/180.0);
    
    // v_d and alt are set to zero for 2d navi framework
    double v_d = 0;
    double alt = 0;
    
    ins_master_struc->timestamp_GPS = GPS_raw_data[0];
    ins_master_struc->hdop = GPS_raw_data[5];
    
    ins_master_struc->processed_data[0] = GPS_raw_data[1];
    ins_master_struc->processed_data[1] = GPS_raw_data[2];
    ins_master_struc->processed_data[2] = alt;
    ins_master_struc->processed_data[3] = v_n;
    ins_master_struc->processed_data[4] = v_e;
    ins_master_struc->processed_data[5] = v_d;
    
}

