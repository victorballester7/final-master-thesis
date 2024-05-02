#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "pointvortices_field.h"
using namespace std;

void saveData(int n, double* X, string filename_output, void* param) {
  ofstream file_output;
  pointvortices_params* prm = (pointvortices_params*)param;
  static uint plot_count = 0;
  // add the plot count to the filename always with 4 digits
  string filename = filename_output + "." + string(4 - to_string(plot_count).length(), '0') + to_string(plot_count) + ".txt";
  file_output.open(filename);
  for (int i = 0; i < n; i++) {
    file_output << X[prm->dim * i] << " " << X[prm->dim * i + 1] << " " << (X[prm->dim * i + 2] > 0) << endl;
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

// computes the square of the velocity of the fluid at the position (x,y) due to the point vortices
double uv2(int n, double* X, double x, double y, void* param) {
  double u = 0;
  double v = 0;
  double xx, yy, r2;
  pointvortices_params* prm = (pointvortices_params*)param;

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
void EnergyProf(int n, double* X, double R_max, string filename_output, void* param) {
  ofstream file_output;
  static int plot_count = 0;
  string filename = filename_output + "." + string(4 - to_string(plot_count).length(), '0') + to_string(plot_count) + ".txt";
  file_output.open(filename);

  double r, E, aux;

  double N_R = 100;      // number of points in the radial direction
  double N_theta = 100;  // number of points in the angular direction

  // for each radius we compute the avergae energy of the system in the ring of radius r
  // E = 1/N_theta sum_{j=1}^{N_theta} (u_j^2 + v_j^2) / 2
  for (int i = 1; i <= N_R; i++) {
    r = R_max * i / N_R;
    E = 0;
    for (int j = 0; j < N_theta; j++) {
      aux = 2 * M_PI * j / N_theta;
      E += uv2(n, X, r * cos(aux), r * sin(aux), param);
    }
    E /= (2 * N_theta);
    file_output << r << " " << E << endl;
  }

  file_output.close();
  plot_count++;
}