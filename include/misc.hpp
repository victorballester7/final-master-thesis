#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void saveData(int n, double* X, double* C, double t, int plot_count, string filename_output) {
  ofstream file_output;
  // add the plot count to the filename always with 4 digits
  string filename = filename_output + "_" + string(4 - to_string(plot_count).length(), '0') + to_string(plot_count) + ".txt";
  file_output.open(filename);
  for (int i = 0; i < n; i++) {
    file_output << X[2 * i] << " " << X[2 * i + 1] << " " << (C[i] > 0) << endl;
  }
  file_output.close();
}

pair<double, double> standard_normal() {
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;

  double z0 = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
  double z1 = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);

  return {z0, z1};
}