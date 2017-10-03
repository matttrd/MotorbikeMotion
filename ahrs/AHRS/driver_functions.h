//
//  driver_functions.h
//  GNSS_INS
//
#ifndef driver_functions_h
#define driver_functions_h

#include <stdio.h>

#include "../../EKF/EKFstruct.h"
#include "../../EKF/EKF.h"
#include "/ekf_functions/INS_matrix_fun.h"
#include "sensor_fusion.h"

#define KNOTS_TO_METERS  0.51444444444444444444444444444444444
#define INIT_STOP_COUNTER 99


void updateGNSS_INS(precision *measurements, EKF *ekf, precision dt,GNSS_INS_t* ins_struct);
void preprocess_GPS_data(double *GPS_raw_data,GNSS_INS_t*);

#endif /* driver_functions_h */
