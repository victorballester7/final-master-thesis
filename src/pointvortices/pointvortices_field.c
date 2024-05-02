#include "../../include/pointvortices_field.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int pointvortices_field(int n, double t, double x[], double f[], void *param) {
  // the field is ordered in x[] and f[] as x1, y1, C1, x2, y2, C2, ..., xn, yn, Cn,
  // where the index i is the particle i-th. Ci are the circulations of the particles
  pointvortices_params *prm = (pointvortices_params *)param;
  // positions
  double twopi = 6.283185307179586;
  double xx, yy, r2;
  int n_bodies = n / prm->dim;
  for (int i = 0; i < n_bodies; i++) {
    f[prm->dim * i] = 0;
    f[prm->dim * i + 1] = 0;
    f[prm->dim * i + 2] = -prm->alpha * x[prm->dim * i + 2];
    for (int k = 0; k < n_bodies; k++) {
      if (k == i) continue;
      xx = x[prm->dim * i] - x[prm->dim * k];
      yy = x[prm->dim * i + 1] - x[prm->dim * k + 1];
      r2 = xx * xx + yy * yy + prm->EPS * prm->EPS;
      f[prm->dim * i] += -x[prm->dim * k + 2] * yy / r2;
      f[prm->dim * i + 1] += x[prm->dim * k + 2] * xx / r2;
    }
    f[prm->dim * i] /= twopi;
    f[prm->dim * i + 1] /= twopi;
  }

  return 0;
}
