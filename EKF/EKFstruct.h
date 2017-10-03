//
//  EKFstruct.h
//  EKF

//#define FIX_R	1;
#define FIX_Q	1;

#define DOUBLE double
#define FLOAT float
#define PRECISION DOUBLE

#ifndef EKFstruct_h
#define EKFstruct_h

typedef PRECISION precision;


typedef struct ekf {
	// state
	int n;	// state dimension
	int m;	// number of measures
	int p;	// number of inputs
	precision *x;	// state
	precision *z;	// observations
	precision *u;
	//precision *u_old;
	precision *y;	// estimated outputs
	
	// functions
	void (*compute_f_F)(struct ekf *, precision dt);
	precision *f;
	// h pointer
	void (*compute_h_H)(struct ekf *);
	precision *h;
	
	// matrices
	precision *F;
	precision *H;
	precision *Q;
	precision *R;
	precision *P;
	
	//update Q and R
	void (*update_Q)(precision *output, precision dt);
	void (*update_R)(struct ekf*);
	
	// additional parameters
	void *parameters;
} EKF;

#endif /* EKFstruct_h */
