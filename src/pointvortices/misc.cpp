#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>

#include "../../include/misc.hpp"
#include "../../include/pointvortices_field.h"
using namespace std;

const uint num_digits_files = 5;

void readData(int *n, double **X, int status, void *param) {
  ifstream file_input;
  string filename_input =
      "data/pointvortices/disk/positions/positions." +
      string(num_digits_files - to_string(status).length(), '0') +
      to_string(status) + ".txt";
  pointvortices_params *prm = (pointvortices_params *)param;
  // file is of the form
  // x1 y1 c1
  // x2 y2 c2
  // ...
  // xn yn cn

  // count the number of lines in the file
  file_input.open(filename_input);
  if (!file_input.is_open()) {
    cout << "Error opening file " << filename_input << endl;
    exit(1);
  }
  *n = 0;
  string line;
  while (getline(file_input, line)) {
    (*n)++;
  }
  file_input.close();

  // allocate memory
  *X = (double *)malloc(prm->dim * (*n) * sizeof(double));
  if (!*X) {
    cout << "Error allocating memory" << endl;
    exit(1);
  }

  // read the data
  file_input.open(filename_input);
  if (!file_input.is_open()) {
    cout << "Error opening file " << filename_input << endl;
    exit(1);
  }
  for (int i = 0; i < *n; i++) {
    file_input >> (*X)[prm->dim * i] >> (*X)[prm->dim * i + 1] >>
        (*X)[prm->dim * i + 2];
  }
  file_input.close();
}

void saveData(int n, double *X, string filename_output, string filename_status,
              double t, void *param, uint plot_count) {
  ofstream file_output;
  pointvortices_params *prm = (pointvortices_params *)param;
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
    file_output << X[prm->dim * i] << " " << X[prm->dim * i + 1] << " "
                << X[prm->dim * i + 2] << endl;
  }
  file_output.close();

  // save the status
  file_output.open(filename_status);
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename_status << endl;
    exit(1);
  }
  file_output << plot_count << endl;
  file_output << t << endl;
  file_output.close();
}

void EnergyTimes(int count_energy, double t, string filename_output) {
  ofstream file_output;
  file_output.open(filename_output, ios::app); // append to the file
  if (!file_output.is_open()) {
    cout << "Error opening file " << filename_output << endl;
    exit(1);
  }
  string status =
      string(num_digits_files - to_string(count_energy).length(), '0') +
      to_string(count_energy);
  file_output << status << " " << t << endl;
  file_output.close();
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

void removeBody(double **X, int j, int *n, void *param) {
  pointvortices_params *prm = (pointvortices_params *)param;
  // remove the j-th body
  // we do this by moving the last body to the j-th position
  // and reducing the number of bodies by one
  for (int k = 0; k < prm->dim; k++) {
    (*X)[prm->dim * j + k] = (*X)[prm->dim * (*n - 1) + k];
  }
  (*n)--;
  double *X_tmp = (double *)realloc(*X, prm->dim * (*n) * sizeof(double));
  if (X_tmp) {
    *X = X_tmp;
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
  (*n)++;
  double *X_tmp = (double *)realloc(*X, prm->dim * (*n) * sizeof(double));
  bool *isInner_tmp = (bool *)realloc(*isInner, *n * sizeof(bool));
  r = sqrt(abs((double)rand() / RAND_MAX)) * R_in;
  if (X_tmp && isInner_tmp) {
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

  double *X_tmp =
      (double *)realloc(*X, prm->dim * (*n + n_extra) * sizeof(double));
  if (X_tmp) {
    *X = X_tmp;
  } else {
    cout << "Error reallocating memory (+)" << endl;
    exit(1);
  }

  for (int i = 0; i < n_extra; i++) {
    if (i % 2 == 0) { // n_extra is even
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

void Energy(int n, double t, double *X, string filename_output, void *param,
            uint plot_count) {
  ofstream file_output;
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
}

// compute the energy of the system in the ring of radius r
void EnergyProf(int n, double *X, int N_R, double *E, double dr,
                string filename_output, void *param, uint plot_count,
                bool save_data) {
  ofstream file_output;
  string filename =
      filename_output + "." +
      string(num_digits_files - to_string(plot_count).length(), '0') +
      to_string(plot_count) + ".txt";

  if (save_data) {
    file_output.open(filename);
    if (!file_output.is_open()) {
      cout << "Error opening file " << filename << endl;
      exit(1);
    }
  }

  double r, aux;

  double N_theta = 1000; // number of points in the angular direction

  // for each radius we compute the avergae energy of the system in the ring
  // of radius r E = 1/N_theta sum_{j=1}^{N_theta} (u_j^2 + v_j^2) / 2
  for (int i = 1; i <= N_R; i++) {
    r = dr * i;
    E[i - 1] = 0;
    for (int j = 0; j < N_theta; j++) {
      aux = 2 * M_PI * j / N_theta;
      E[i - 1] += uv2(n, X, r * cos(aux), r * sin(aux), param);
    }
    E[i - 1] /= (2 * N_theta);
    if (save_data)
      file_output << r << " " << E[i - 1] << endl;
  }

  if (save_data) {
    file_output.close();
  }
}

// compute the energy flux of the system in disks of radius r
void EnergyFlux(double dt, int N_R, double *E0, double *E1, double dr,
                string filename_output, uint plot_count) {
  ofstream file_output;
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
    r = dr * i;
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
}

void NumVortices(int n, double *X, int N_R, double dr, string filename_output,
                 void *param, uint plot_count) {
  ofstream file_output;
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
  int aux;
  for (int i = 0; i < n; i++) {
    r = sqrt(X[prm->dim * i] * X[prm->dim * i] +
             X[prm->dim * i + 1] * X[prm->dim * i + 1]);
    aux = (int)(r / dr);
    if (aux < N_R)
      count[aux]++;
  }
  for (int i = 1; i <= N_R; i++) {
    r = dr * i;
    file_output << r << " " << count[i - 1] << " " << n << endl;
  }

  file_output.close();
  delete[] count;
}
