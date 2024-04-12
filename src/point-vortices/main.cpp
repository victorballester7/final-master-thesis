#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../include/misc.hpp"
#include "../../include/pointvortices_field.h"
#include "../../include/rk78.h"

#define TOL 1e-6  // tolerance for the RK78 method
#define RK78 1
using namespace std;

int main(int argc, char* argv[]) {
  const string filename_input = "config/pointvortices/input.txt";  // name of the input file to read the parameters
  const string filename_output = "data/pointvortices/positions";   // name of the output file to write the positions of the point vortices
  const string space = "    ";                                     // space to print the parameters
  const uint per = 10;                                             // progress percentage interval to print (each count%)
  pointvortices_params prm;                                        // parameters of the system
  prm.dim = 2;                                                     // dimension of the space
  int n;                                                           // number of bodies
  uint count = per;
  double plot_dt;       // frequency to write the animation file
  double R;             // radius of the circle in which the bodies are placed
  double C;             // circulation
  double dt;            // time step
  double T;             // final time
  uint plot_count = 0;  // counter to write the animation file

  // ------------- File input setup ----------------
  ifstream file_input;
  ofstream file_output;
  file_input.open(filename_input);
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> n;
    file_input >> tmp >> C;
    file_input >> tmp >> R;
    file_input >> tmp >> dt;
    file_input >> tmp >> T;
    file_input >> tmp >> plot_dt;
  }
  file_input.close();

  // print parameters
  // ------------- Print plot setup -----------------
  cout << "Number of votices:  " << space << n << endl;
  cout << "dt:                 " << space << dt << endl;
  cout << "Tfinal:             " << space << T << endl;
  cout << "Initial radius:     " << space << R << endl;
  cout << "Circulation:        " << space << C << endl;

  const int N = prm.dim * n;  // dimension of the field
  prm.EPS = 1e-3;             // softening parameter
  // -------------------------------------------

  double t = 0.0;
  uint numSteps = 0;
  double fraction_completed = T / 100.;  // fraction of the integration time to print

  // allocate memory
  double* X = new double[N];  // vector of positions (x1,y1,x2,y2,...,xn,yn)
  prm.C = new double[n];      // circulations

  // set circulations
  double eps = C / 100;  // small perturbation to avoid the same circulations
  pair<double, double> z;
  for (int i = 0; i < n; i += 2) {
    z = standard_normal();
    prm.C[i] = C + eps * z.first;
    prm.C[i + 1] = -C + eps * z.second;
  }

  // set initial conditions (random in the circle of radius R)
  double aux;
  for (int i = 0; i < N; i += 2) {
    aux = (double)rand() / RAND_MAX;
    X[i] = R * cos(2 * M_PI * aux);
    X[i + 1] = R * sin(2 * M_PI * aux);
  }

  // save initial conditions
  saveData(n, X, prm.C, t, plot_count, filename_output);

  plot_count++;

  numSteps++;

  double hmin = dt / 100;
  double hmax = dt * 100;

  while (t < T - TOL) {
    // simulate

#ifdef RK78
    if (rk78(&t, X, &dt, hmin, hmax, TOL, N, pointvortices_field, &prm) != 0) {
      cout << "Error in the integration" << endl;
      return 1;
    }
#else
    if (rk4(&t, X, &dt, N, pointvortices_field, &prm) != 0) {
      cout << "Error in the integration" << endl;
      return 1;
    }
#endif
    // save data
    if (t > plot_count * plot_dt - TOL) {
      saveData(n, X, prm.C, t, plot_count, filename_output);
      plot_count++;
    }

    // print progress
    if (t > fraction_completed * count - TOL) {
      cout << count << "%"
           << " dt: " << dt << endl;
      count += per;
    }

    numSteps++;
  }

  return 0;
}
