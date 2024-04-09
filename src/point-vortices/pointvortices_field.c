#include "../../include/pointvortices_field.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../include/rk78.h"

int pointvortices_field(int n, double t, double x[], double f[], void *param) {
  // the field is ordered in x[] and f[] as x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, vx2, vy2, vz2, ..., where the index i is the particle i-th. If the field is 2D, then the same order is used but without the z components.
  pointvortices_params *prm = (pointvortices_params *)param;
  // positions
  double twopi = 6.283185307179586;
  double xx, yy, r2;
  int n_bodies = n / prm->dim;
  for (int i = 0; i < n_bodies; i++) {
    f[2 * i] = 0;
    f[2 * i + 1] = 0;
    for (int k = 0; k < n_bodies; k++) {
      if (k == i) continue;
      xx = x[2 * i] - x[2 * k];
      yy = x[2 * i + 1] - x[2 * k + 1];
      r2 = xx * xx + yy * yy + prm->EPS * prm->EPS;
      f[2 * i] += -prm->C[k] * yy / r2;
      f[2 * i + 1] += prm->C[k] * xx / r2;
    }
    f[2 * i] /= twopi;
    f[2 * i + 1] /= twopi;
  }

  return 0;
}
