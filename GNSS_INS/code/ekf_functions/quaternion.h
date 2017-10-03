//
//  quaternion.h

#ifndef quaternion_h
#define quaternion_h

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>



typedef struct quaternion quaternion_t;

struct quaternion
{
    double q[4];
};

quaternion_t *quaternion_new(void);

quaternion_t *quaternion_new_set(double,
                                 double,
                                 double,
                                 double);

void quaternion_copy(quaternion_t*, quaternion_t*);

double quaternionNorm(quaternion_t*);

void quaternConj(quaternion_t*, quaternion_t*);

bool quaternionIsEqual(quaternion_t*, quaternion_t*);

void quaternProd(quaternion_t*, quaternion_t*, quaternion_t*);


void quaternion_print(quaternion_t*);

#endif 
