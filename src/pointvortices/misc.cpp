#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "../../include/misc.hpp"
#include "../../include/pointvortices_field.h"
using namespace std;

const uint num_digits_files = 5;

void saveData(int n, double *X, string filename_output, void *param) {
  ofstream file_output;
  pointvortices_params *prm = (pointvortices_params *)param;
  static uint plot_count = 0;
  // add the plot count to the filename always with 4 digits
  string filename =
      filename_output + "." +
      string(num_digits_files - to_string(plot_count).length(), '0') +
      to_string(plot_count) + ".txt";
  file_output.open(filename);
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename << endl;
    exit(1);
  }
  for (int i = 0; i < n; i++) {
    int j;
    if (X[prm->dim * i] * X[prm->dim * i] +
            X[prm->dim * i + 1] * X[prm->dim * i + 1] <
        0.25 * 0.25) {
      j = 11111;
    } else {
      j = 22222;
    }
    file_output << X[prm->dim * i] << " " << X[prm->dim * i + 1] << " "
                << (X[prm->dim * i + 2] > 0) << " " << j << endl;
  }
  file_output.close();
  plot_count++;
}

pair<double, double> standard_normal() {
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;

  double z0 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
  double z1 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);

  return {z0, z1};
}

void resetBody(double *X, int i, double R, void *param) {
  pointvortices_params *prm = (pointvortices_params *)param;
  double r, aux;
  aux = 2 * M_PI * (double)rand() / RAND_MAX;
  r = sqrt(abs((double)rand() / RAND_MAX)) * R;
  X[prm->dim * i] = r * cos(aux);
  X[prm->dim * i + 1] = r * sin(aux);
}

void removeBody(double **X, bool **isInner, int j, int *n, void *param) {
  pointvortices_params *prm = (pointvortices_params *)param;
  // remove the j-th body
  // we do this by moving the last body to the j-th position
  // and reducing the number of bodies by one
  for (int k = 0; k < prm->dim; k++) {
    (*X)[prm->dim * j + k] = (*X)[prm->dim * (*n - 1) + k];
  }
  (*isInner)[j] = (*isInner)[*n - 1];
  (*n)--;
  double *X_tmp;
  bool *isInner_tmp;
  if ((X_tmp = (double *)realloc(*X, prm->dim * (*n) * sizeof(double))) &&
      (isInner_tmp = (bool *)realloc(*isInner, *n * sizeof(bool)))) {
    *X = X_tmp;
    *isInner = isInner_tmp;
  } else {
    cout << "Error reallocating memory (-)" << endl;
    exit(1);
  }
}

void addBody(double **X, bool **isInner, int *n, double R_in, double C,
             void *param) {
  pointvortices_params *prm = (pointvortices_params *)param;
  double r, aux;
  aux = 2 * M_PI * (double)rand() / RAND_MAX;
  double *X_tmp;
  bool *isInner_tmp;
  r = sqrt(abs((double)rand() / RAND_MAX)) * R_in;
  (*n)++;
  if ((X_tmp = (double *)realloc(*X, prm->dim * (*n) * sizeof(double))) &&
      (isInner_tmp = (bool *)realloc(*isInner, *n * sizeof(bool)))) {
    *X = X_tmp;
    *isInner = isInner_tmp;
  } else {
    cout << "Error reallocating memory (+)" << endl;
    exit(1);
  }

  (*X)[prm->dim * (*n - 1)] = r * cos(aux);
  (*X)[prm->dim * (*n - 1) + 1] = r * sin(aux);
  (*X)[prm->dim * (*n - 1) + 2] = C;
  (*isInner)[*n - 1] = true;
}

void addVortices(double **X, int *n, int n_extra, double R_in, double eps,
                 void *param) {
  pair<double, double> z;
  pointvortices_params *prm = (pointvortices_params *)param;
  double r, aux;

  double *X_tmp;
  if ((X_tmp =
           (double *)realloc(*X, prm->dim * (*n + n_extra) * sizeof(double)))) {
    *X = X_tmp;
  } else {
    cout << "Error reallocating memory (+)" << endl;
    exit(1);
  }

  for (int i = 0; i < n_extra; i++) {
    if (i % 2 == 0) {
      z = standard_normal();
      (*X)[prm->dim * (*n + i) + 2] = eps * z.first;
      (*X)[prm->dim * (*n + i + 1) + 2] = -eps * z.first;
    }
    aux = 2 * M_PI * (double)rand() / RAND_MAX;
    r = R_in * sqrt(abs((double)rand() /
                        RAND_MAX)); // random radius. We do not do it uniformly
                                    // to avoid having more points per area in
                                    // the smaller subradii. With the sqrt we
                                    // have more points in the outer layers.
    (*X)[prm->dim * (*n + i)] = r * cos(aux);     // x position
    (*X)[prm->dim * (*n + i) + 1] = r * sin(aux); // y position
  }
  *n += n_extra;
}

// computes the square of the velocity of the fluid at the position (x,y)
// due to the point vortices
double uv2(int n, double *X, double x, double y, void *param) {
  double u = 0;
  double v = 0;
  double xx, yy, r2;
  pointvortices_params *prm = (pointvortices_params *)param;

  for (int i = 0; i < n; i++) {
    xx = x - X[prm->dim * i];
    yy = y - X[prm->dim * i + 1];
    r2 = xx * xx + yy * yy + prm->EPS * prm->EPS;
    u += -X[prm->dim * i + 2] * yy / r2;
    v += X[prm->dim * i + 2] * xx / r2;
  }
  u /= 2 * M_PI;
  v /= 2 * M_PI;
  return u * u + v * v;
}

void Energy(int n, double t, double *X, string filename_output, void *param) {
  ofstream file_output;
  static int plot_count = 0;
  file_output.open(filename_output, ios::app); // append to the file
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename_output << endl;
    exit(1);
  }

  pointvortices_params *prm = (pointvortices_params *)param;
  double E = 0;
  double xx, yy, r2;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      xx = X[prm->dim * i] - X[prm->dim * j];
      yy = X[prm->dim * i + 1] - X[prm->dim * j + 1];
      r2 = xx * xx + yy * yy + prm->EPS * prm->EPS;
      E += X[prm->dim * i + 2] * X[prm->dim * j + 2] * log(r2);
    }
  }
  E /= 2 * M_PI;

  file_output << t << " " << E << endl;
  file_output.close();
  plot_count++;
}

// compute the energy of the system in the ring of radius r
void EnergyProf(int n, double *X, int N_R, double *E, double R_max,
                string filename_output, void *param) {
  ofstream file_output;
  static int plot_count = 0;
  string filename =
      filename_output + "." +
      string(num_digits_files - to_string(plot_count).length(), '0') +
      to_string(plot_count) + ".txt";

  file_output.open(filename);
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename << endl;
    exit(1);
  }

  double r, aux;

  double N_theta = 1000; // number of points in the angular direction

  // for each radius we compute the avergae energy of the system in the ring
  // of radius r E = 1/N_theta sum_{j=1}^{N_theta} (u_j^2 + v_j^2) / 2
  for (int i = 1; i <= N_R; i++) {
    r = R_max * i / N_R;
    E[i - 1] = 0;
    for (int j = 0; j < N_theta; j++) {
      aux = 2 * M_PI * j / N_theta;
      E[i - 1] += uv2(n, X, r * cos(aux), r * sin(aux), param);
    }
    E[i - 1] /= (2 * N_theta);
    file_output << r << " " << E[i - 1] << endl;
  }

  file_output.close();
  plot_count++;
}

// compute the energy flux of the system in disks of radius r
void EnergyFlux(double dt, int N_R, double *E0, double *E1, double R_max,
                string filename_output) {
  ofstream file_output;
  static int plot_count = 0;
  string filename =
      filename_output + "." +
      string(num_digits_files - to_string(plot_count).length(), '0') +
      to_string(plot_count) + ".txt";

  file_output.open(filename);
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename << endl;
    exit(1);
  }

  double flux, r;
  double E1_total, E0_total;
  for (int i = 1; i <= N_R; i++) {
    r = R_max * i / N_R;
    E1_total = 0;
    E0_total = 0;
    // compute the total energy inside the circle of radius r
    for (int j = 0; j < i; j++) {
      E1_total += E1[j];
      E0_total += E0[j];
    }
    flux = (E1_total - E0_total) / dt;
    file_output << r << " " << flux << endl;
  }
  file_output.close();
  plot_count++;
}

void NumVortices(int n, double *X, int N_R, double R_max,
                 string filename_output, void *param) {
  ofstream file_output;
  static int plot_count = 0;
  string filename =
      filename_output + "." +
      string(num_digits_files - to_string(plot_count).length(), '0') +
      to_string(plot_count) + ".txt";
  file_output.open(filename);
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename << endl;
    exit(1);
  }

  pointvortices_params *prm = (pointvortices_params *)param;

  double *count = new double[N_R];
  for (int i = 0; i < N_R; i++)
    count[i] = 0;
  double r;
  for (int i = 0; i < n; i++) {
    r = sqrt(X[prm->dim * i] * X[prm->dim * i] +
             X[prm->dim * i + 1] * X[prm->dim * i + 1]);
    count[(int)(r / R_max * N_R)]++;
  }
  for (int i = 1; i <= N_R; i++) {
    file_output << R_max * i / N_R << " " << count[i - 1] << " " << n << endl;
  }

  file_output.close();
  plot_count++;
  delete[] count;
}
