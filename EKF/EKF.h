//
//  EKF.h
//  EKF

#ifndef EKF_h
#define EKF_h

#include <stdio.h>
#include "EKFstruct.h"

//----------------------------------------- LINEAR ALGEBRA UTILITIES --------------------------------------
static inline void transpose(precision *output, int m, int n, precision *H);

static inline void symmetrize(int n, precision *mat);

//---------------------------------------------------------------------------------------------------------

EKF *init_ekf(int n, int m, int p, precision *init_state, precision *init_Q, precision *init_R, precision *init_P);

void free_ekf(EKF *ekf);

void ekf_prediction(EKF *ekf, precision dt);

void ekf_update(EKF *ekf);

#endif /* EKF_h */

