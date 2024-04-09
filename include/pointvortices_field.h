#ifndef N_BODY_FIELD_H
#define N_BODY_FIELD_H

#ifdef __cplusplus
extern "C" {
#endif
typedef struct pointvortices_params {
  int dim;     // dimension of the space (2 or 3)
  double *C;   // circulations
  double EPS;  // softening parameter
} pointvortices_params;

// -----------------------------------------------------
// pointvortices_field
// -----------------------------------------------------
// Purpose:
// 	Compute the n-body field at a given point
//
// Parameters:
// 	n: dimension of the field = 2 * 2 * n_bodies (in R^2) or 3 * 2 * n_bodies (in R^3)
// 	t: time at the integration preces (not used in the equations, but necessary to pass to the RK78)
// 	x: point at which the field is evaluated
// 	f: output field
// 	param: pointer to the parameters of the system (masses, gravitational constant)
//
// Returns:
// 	0 if everything went fine
//  1 if otherwise
// -----------------------------------------------------
int pointvortices_field(int n, double t, double x[], double f[], void *param);

#ifdef __cplusplus
}
#endif

#endif
