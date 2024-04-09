#ifndef RK78_H
#define RK78_H

#ifdef __cplusplus
extern "C" {
#endif
#define SGN(x) ((x) >= 0 ? 1 : -1)  // sign function

// ----------------------------------------------
// rk78 coefficients
// ----------------------------------------------
// const double alfa[13] = {
//     0., 2. / 27., 1. / 9., 1. / 6., 5. / 12., .5, 5. / 6., 1. / 6., 2. / 3., 1. / 3., 1., 0., 1.};
// const double beta[78] = {
//     2. / 27., 1. / 36., 1. / 12., 1. / 24., 0., 1. / 8., 5. / 12., 0., -25. / 16., 25. / 16., .05, 0., 0., .25, .2, -25. / 108., 0., 0., 125. / 108., -65. / 27., 2. * (125. / 108.), 31. / 300., 0., 0., 0., 61. / 225., -2. / 9., 13. / 900., 2., 0., 0., -53. / 6., 704. / 45., -107. / 9., 67. / 90., 3., -91. / 108., 0., 0., 23. / 108., -976. / 135., 311. / 54., -19. / 60., 17. / 6., -1. / 12., 2383. / 4100., 0., 0., -341. / 164., 4496. / 1025., -301. / 82., 2133. / 4100., 45. / 82., 45. / 164., 18. / 41., 3. / 205., 0., 0., 0., 0., -6. / 41., -3. / 205., -3. / 41., 3. / 41., 6. / 41., 0., -1777. / 4100., 0., 0., -341. / 164., 4496. / 1025., -289. / 82., 2193. / 4100., 51. / 82., 33. / 164., 12. / 41., 0., 1.};
// const double c[11] = {41. / 840., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 41. / 840.};
// const double cp[13] = {0., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 0., 41. / 840., 41. / 840.};

// ----------------------------------------------
// rk78
// ----------------------------------------------
// Purpose:
//    Runge-Kutta-Fehlberg 7(8) integrator
//
// Parameters:
//    t: pointer to the current time (initially pointing to the time of the initial conditions and at the end of the function pointing to the time at t + h)
//    x: pointer to the current state vector (initially pointing to the initial conditions and at the end of the function pointing to the state vector at t + h)
//    h: pointer to the wished step size (initially the step at which we want to integrate the system and at the end of the function pointing to the real step size of the next step)
//    hmin: minimum step size (if the step size is < 0, then hmin < 0)
//    hmax: maximum step size (if the step size is < 0, then hmax < 0)
//    tol: tolerance, which corresponds to the (theoretical) error between the real solution and the next step of the solution
//    n: dimension of the field
//    field: pointer to the function that calculates the field
//    param: pointer to the parameters of the field function
//
// Return value:
//    0: success
//    otherwise: error
// ----------------------------------------------
int rk78(double *t, double x[], double *h, double hmin, double hmax, double tol, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param);

int rk4(double *t, double x[], double *h, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param);
// ----------------------------------------------
// flow
// ----------------------------------------------
// Purpose:
//    Integrates the field from t to t + T
//
// Parameters:
//    t: pointer to the current time (initially pointing to the time of the initial conditions and at the end of the function pointing to the time at t + T)
//    x: pointer to the current state vector (initially pointing to the initial conditions and at the end of the function pointing to the state vector at t + T)
//    h: pointer to the wished step size (initially the solution time of the next step will as as close as possible to *t + h, where *t is the inital time)
//    T: integration time
//    hmin: minimum step size
//    hmax: maximum step size
//    tol: tolerance, which corresponds to the (theoretical) error between the real solution and the next step of the solution
//    maxNumSteps: maximum number of steps
//    n: dimension of the field
//    field: pointer to the function that calculates the field
//    param: pointer to the parameters of the field function
//
// Return value:
//    0: success
//    otherwise: error
// ----------------------------------------------
int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param);

#ifdef __cplusplus
}
#endif

#endif  // RK78_H
