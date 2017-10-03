#ifndef initialization_h
#define initialization_h

#include "../../EKF/EKFstruct.h"

// Initial stage functions
void initialization_stage_estimates(EKF *ekf, precision *measurements, int init_counter);

void GNSS_INS_initialization(EKF *ekf,GNSS_INS_t *ins_struct);
#endif
