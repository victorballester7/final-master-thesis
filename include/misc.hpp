#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void saveData(int n, double* X, double* C, string filename_output) {
  ofstream file_output;
  static uint plot_count = 0;
  // add the plot count to the filename always with 4 digits
  string filename = filename_output + "." + string(4 - to_string(plot_count).length(), '0') + to_string(plot_count) + ".txt";
  file_output.open(filename);
  for (int i = 0; i < n; i++) {
    file_output << X[2 * i] << " " << X[2 * i + 1] << " " << (C[i] > 0) << endl;
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

// compute the energy the kinetic energy of the body i influenced by the rest of the bodies
double energy(int n, double* X, double* C, int body) {
  double E = 0, dist;
  for (int j = 0; j < n; j++) {
    if (j == body) continue;
    dist = (X[2 * body] - X[2 * j]) * (X[2 * body] - X[2 * j]) + (X[2 * body + 1] - X[2 * j + 1]) * (X[2 * body + 1] - X[2 * j + 1]);
    E += C[j] * log(dist);
  }
  E *= -C[body] / (4 * M_PI);
  return E;
}

// compute the energy of the system in
void EnergyProf(int n, double* X, double* C, double R_max, string filename_output) {
  ofstream file_output;
  static int plot_count = 0;
  string filename = filename_output + "." + string(4 - to_string(plot_count).length(), '0') + to_string(plot_count) + ".txt";
  file_output.open(filename);

  double r, E, aux;

  double N_R = 100;  // number of points in the radial direction
  uint count;

  for (int i = 0; i < N_R; i++) {
    r = R_max * i / N_R;
    E = 0;
    count = 0;
    for (int j = 0; j < n; j++) {
      aux = (X[2 * j] * X[2 * j] + X[2 * j + 1] * X[2 * j + 1]);
      if (aux < r * r) {
        E += energy(n, X, C, j);
        count++;
      }
    }
    if (count > 0) {
      E /= count;
      file_output << r << " " << E << endl;
    }
  }

  file_output.close();
  plot_count++;
}