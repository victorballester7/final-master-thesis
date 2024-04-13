#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const double c[13] = {
    0., 2. / 27., 1. / 9., 1. / 6., 5. / 12., .5, 5. / 6., 1. / 6., 2. / 3., 1. / 3., 1., 0., 1.};
const double A[78] = {
    2. / 27., 1. / 36., 1. / 12., 1. / 24., 0., 1. / 8., 5. / 12., 0., -25. / 16., 25. / 16., .05, 0., 0., .25, .2, -25. / 108., 0., 0., 125. / 108., -65. / 27., 2. * (125. / 108.), 31. / 300., 0., 0., 0., 61. / 225., -2. / 9., 13. / 900., 2., 0., 0., -53. / 6., 704. / 45., -107. / 9., 67. / 90., 3., -91. / 108., 0., 0., 23. / 108., -976. / 135., 311. / 54., -19. / 60., 17. / 6., -1. / 12., 2383. / 4100., 0., 0., -341. / 164., 4496. / 1025., -301. / 82., 2133. / 4100., 45. / 82., 45. / 164., 18. / 41., 3. / 205., 0., 0., 0., 0., -6. / 41., -3. / 205., -3. / 41., 3. / 41., 6. / 41., 0., -1777. / 4100., 0., 0., -341. / 164., 4496. / 1025., -289. / 82., 2193. / 4100., 51. / 82., 33. / 164., 12. / 41., 0., 1.};
const double brk7[11] = {41. / 840., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 41. / 840.};
const double brk8[13] = {0., 0., 0., 0., 0., 34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280., 0., 41. / 840., 41. / 840.};

int rk78(double *t, double x[], double *h,
         double hmin, double hmax, double tol,
         int n, int (*field)(int n, double t, double x[], double f[], void *param),
         void *param) {
  // Local variables
  int return_value;  // to store the return value of field
  int indexA;        // to store the index of A
  double tt, dif, RK8norm, eps, temp;
  double K[13 * n];  // matrix vector containing all the stages
  double RK7[n];     // vector containing the RK7 step
  double RK8[n];     // vector containing the RK8 step

  do {  // do a step for both the RK7 and RK8 methods
    indexA = 0;
    // loop over stages (RK7 has 11 stages and RK8 has 13 stages)
    for (int i = 0; i < 13; i++) {
      memcpy(RK7, x, n * sizeof(double));  // we use RK7 as a temporary vector (in order to avoid using more memory)
      tt = (*t) + c[i] * (*h);             // time at stage i
      for (int j = 0; j < i; j++, indexA++) {
        temp = A[indexA] * (*h);     // temp = step size * A[i][j]
        for (int m = 0; m < n; m++)  // sum over all the components of the contribution of stage j
          RK7[m] += temp * K[j * n + m];
      }
      return_value = field(n, tt, RK7, K + i * n, param);
      if (return_value) return return_value;  // Error
    }
    RK8norm = dif = 0;
    memcpy(RK8, x, n * sizeof(double));
    memcpy(RK7, x, n * sizeof(double));

    // we loop over the coordinates and then over the stages to save time and compute d and dd at the same time.
    for (int m = 0; m < n; m++) {     // loop over the coordinates
      for (int i = 0; i < 11; i++) {  // loop over the first 11 stages that both RK7 and RK8 'share'
        temp = (*h) * K[i * n + m];
        RK7[m] += temp * brk7[i];
        RK8[m] += temp * brk8[i];
      }
      // add the contribution of the last two stages (on RK8)
      RK8[m] += (*h) * (brk8[11] * K[11 * n + m] + brk8[12] * K[12 * n + m]);
      dif += fabs(RK8[m] - RK7[m]);
      RK8norm += fabs(RK8[m]);  // L1 norm of RK8 used to relativize the tolerance
    }
    dif /= n;                          // average of the differences between RK7 and RK8
    eps = tol * (1. + .01 * RK8norm);  // updated tolerance that accounts for large values of the field (in that case the tolerance is increased).
    // if the error is small enough (less than our 'artifical' tolerance) or the step size is smaller than the minimum step size we break the loop
    if (dif < eps || fabs(*h) <= fabs(hmin)) break;
    (*h) *= .9 * pow(eps / dif, .125);  // Fehlberg correction. Note that 0.125 = 1/(1 + 7). See Stoer 3rd edition 7.2.5.13. We take the conservative correction coefficient as 0.9.
    if (fabs(*h) < fabs(hmin)) *h = hmin;
  } while (1);
  (*t) += (*h);                        // we update the time
  memcpy(x, RK8, n * sizeof(double));  // we update the vector x with the RK8 step
  // If the error is super small, we don't want the next step to be huge. So we bound the error from below by eps / 256.
  if (dif < eps / 256) dif = eps / 256;

  (*h) *= .9 * pow(eps / dif, .125);  // Fehlberg correction. Note that 0.125 = 1/(1 + 7). See Stoer 3rd edition 7.2.5.13. We take the conservative correction coefficient as 0.9.
  // we impose the limits on the step size
  if (fabs(*h) < fabs(hmin))
    *h = hmin;
  else if (fabs(*h) > fabs(hmax))
    *h = hmax;
  return 0;
}

int rk4(double *t, double x[], double *h, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param) {
  double k1[n], k2[n], k3[n], k4[n];
  double t0 = *t;

  // stage 1
  if (field(n, t0, x, k1, param)) return 1;

  // stage 2
  for (int i = 0; i < n; i++) {
    x[i] += k1[i] * (*h) / 2;
  }
  *t = t0 + *h / 2;
  if (field(n, *t, x, k2, param)) return 1;

  // stage 3
  for (int i = 0; i < n; i++) {
    x[i] += -k1[i] * (*h) / 2 + k2[i] * (*h) / 2;
  }
  if (field(n, *t, x, k3, param)) return 1;

  // stage 4
  for (int i = 0; i < n; i++) {
    x[i] += -k2[i] * (*h) / 2 + k3[i] * (*h);
  }
  *t = t0 + *h;
  if (field(n, *t, x, k4, param)) return 1;

  for (int i = 0; i < n; i++) {
    x[i] += -k3[i] * (*h) + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * (*h) / 6;
  }
  *t = t0 + *h;
  return 0;
}

int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param) {
  double t0 = *t;
  int count = 0;

  if (fabs(T) < tol) {
    *t = t0 + T;
    return 0;
  }

  while (fabs(*t - t0) < fabs(T) && count < maxNumSteps) {
    if (fabs(*t + *h - t0) > fabs(T)) {
      *h = t0 + T - *t;
      if (rk78(t, x, h, hmin, hmax, tol, n, field, param))
        return 1;
      else
        break;
    }
    if (rk78(t, x, h, hmin, hmax, tol, n, field, param)) return 1;
    count++;
  }
  if (count == maxNumSteps) {
    return 1;
  } else {
    *t = t0 + T;
    return 0;
  }
}
