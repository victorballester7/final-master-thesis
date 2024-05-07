#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "../../include/misc.hpp"
#include "../../include/pointvortices_field.h"
using namespace std;

void saveData(int n, double *X, string filename_output, void *param) {
  ofstream file_output;
  pointvortices_params *prm = (pointvortices_params *)param;
  static uint plot_count = 0;
  // add the plot count to the filename always with 4 digits
  string filename = filename_output + "." +
                    string(4 - to_string(plot_count).length(), '0') +
                    to_string(plot_count) + ".txt";
  file_output.open(filename);
  for (int i = 0; i < n; i++) {
    file_output << X[prm->dim * i] << " " << X[prm->dim * i + 1] << " "
                << (X[prm->dim * i + 2] > 0) << endl;
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
  r = sqrt(R * abs((double)rand() / RAND_MAX));
  X[prm->dim * i] = r * cos(aux);
  X[prm->dim * i + 1] = r * sin(aux);
}

// computes the square of the velocity of the fluid at the position (x,y) due to
// the point vortices
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
// compute the energy of the system in
void EnergyProf(int n, double *X, int N_R, double *E, double R_max,
                string filename_output, void *param) {
  ofstream file_output;
  static int plot_count = 0;
  string filename = filename_output + "." +
                    string(4 - to_string(plot_count).length(), '0') +
                    to_string(plot_count) + ".txt";
  file_output.open(filename);

  double r, aux;

  double N_theta = 1000; // number of points in the angular direction

  // for each radius we compute the avergae energy of the system in the ring of
  // radius r E = 1/N_theta sum_{j=1}^{N_theta} (u_j^2 + v_j^2) / 2
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

void EnergyFlux(double dt, int N_R, double *E0, double *E1, double R_max,
                string filename_output) {
  ofstream file_output;
  static int plot_count = 0;
  string filename = filename_output + "." +
                    string(4 - to_string(plot_count).length(), '0') +
                    to_string(plot_count) + ".txt";
  file_output.open(filename);
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
