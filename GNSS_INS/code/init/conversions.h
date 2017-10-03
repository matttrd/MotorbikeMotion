// conversion functions declarations
#ifndef conversions_h
#define conversions_h

#include <stdlib.h>
#include <math.h>
#include "EKFstruct.h"

void euler_angles_from_acc_and_mag(precision ax, precision ay, precision az, precision mx, precision my, precision mz, precision *ret);

void euler_angles_to_rotation(precision phi, precision theta, precision psi, precision *mat);


void rotation_to_quaternion(precision *R, precision *ret);
#endif