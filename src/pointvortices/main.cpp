#include <cmath>
#include <fstream>
#include <iostream>

#include "../../include/misc.hpp"
#include "../../include/pointvortices_field.h"
#include "../../include/rk78.h"

#define TOL 1e-6  // tolerance for the RK78 method
#define RK78 1
using namespace std;

int main(void) {
  const string filename_input = "config/pointvortices/input.txt";               // name of the input file to read the parameters
  const string filename_output = "data/pointvortices/positions/positions";      // name of the output file to write the positions of the point vortices
  const string filename_output_E = "data/pointvortices/EnergyProf/EnergyProf";  // name of the output file to write the positions of the point vortices
  const string space = "    ";                                                  // space to print the parameters
  pointvortices_params prm;                                                     // parameters of the system
  prm.dim = 3;                                                                  // dimension of the space (2 space dimensions + 1 circulation dimension)
  int n;                                                                        // number of bodies

  double R;          // radius of the circle in which the bodies are placed
  double R_exit;     // radius of the circle in which the bodies are simulated
  double C0;         // initial circulation
  double dt;         // time step
  uint totalSteps;   // total number of steps
  uint outSteps;     // number of steps to write the positions of the bodies
  uint energySteps;  // number of steps to write the energy profile
  uint printSteps;   // number of steps to print the progress

  // ------------- File input setup ----------------
  ifstream file_input;
  ofstream file_output;
  file_input.open(filename_input);
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> totalSteps;
    file_input >> tmp >> outSteps;
    file_input >> tmp >> energySteps;
    file_input >> tmp >> printSteps;
    file_input >> tmp >> n;
    file_input >> tmp >> C0;
    file_input >> tmp >> prm.alpha;
    file_input >> tmp >> R;
    file_input >> tmp >> R_exit;
    file_input >> tmp >> dt;
  }
  file_input.close();

  // print parameters
  // ------------- Print plot setup -----------------
  cout << "Number of steps:    " << space << totalSteps << endl;
  cout << "Number of votices:  " << space << n << endl;
  cout << "Initial radius:     " << space << R << endl;
  cout << "Exit radius:        " << space << R_exit << endl;
  cout << "Circulation:        " << space << C0 << endl;
  cout << "Drag coefficient:   " << space << prm.alpha << endl;

  const int N = prm.dim * n;  // dimension of the field
  prm.EPS = 1e-3;             // softening parameter
  // -------------------------------------------

  double t = 0.0;
  uint numSteps = 0;

  // allocate memory
  double* X = new double[N];  // vector of positions (x1,y1,C1,x2,y2,C2,...,xn,yn,Cn)
  double min_dist, dist;
  int closest;

  // set circulations
  double eps = C0 / 100;  // small perturbation to avoid the same circulations
  pair<double, double> z;
  for (int i = 0; i < n; i += 2) {
    z = standard_normal();  // generates 2 iid standard normal variables using the Box-Muller method
    X[prm.dim * i + 2] = C0 + eps * z.first;
    X[prm.dim * (i + 1) + 2] = -C0 + eps * z.second;
  }

  // seed the random number generator
  srand(time(NULL));

  // set initial conditions (random in the circle of radius R)
  double aux, r;
  for (int i = 0; i < n; i++) {
    aux = 2 * M_PI * (double)rand() / RAND_MAX;
    r = sqrt(R * abs((double)rand() / RAND_MAX));  // random radius. We do not do it uniformly to avoid having more points per area in the smaller subradii. With the sqrt we have more points in the outer layers.
    X[prm.dim * i] = r * cos(aux);                 // x position
    X[prm.dim * i + 1] = r * sin(aux);             // y position
  }

  double hmin = dt / 100;
  double hmax = dt * 100;

  while (numSteps < totalSteps) {
    // save data
    if (numSteps % outSteps == 0) {
      saveData(n, X, filename_output, &prm);
    }

    // plot energy profile
    if (numSteps % energySteps == 0) {
      EnergyProf(n, X, R_exit, filename_output_E, &prm);
    }

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

    // check if a body has exited the domain
    for (int i = 0; i < n; i++) {
      if (X[prm.dim * i] * X[prm.dim * i] + X[prm.dim * i + 1] * X[prm.dim * i + 1] > R_exit * R_exit) {
        // search the closest point to i-th body (advection of 2 bodies in a straight line)
        min_dist = R_exit;  // Radius of the circle
        closest = -1;
        for (int j = 0; j < n; j++) {
          if (i == j) continue;  // do not compare with itself
          dist = (X[prm.dim * i] - X[prm.dim * j]) * (X[prm.dim * i] - X[prm.dim * j]) + (X[prm.dim * i + 1] - X[prm.dim * j + 1]) * (X[prm.dim * i + 1] - X[prm.dim * j + 1]);
          if (dist < min_dist) {
            min_dist = dist;
            closest = j;
          }
        }
        // set again randomly in the circle of radius R the body i
        aux = 2 * M_PI * (double)rand() / RAND_MAX;
        r = sqrt(R * abs((double)rand() / RAND_MAX));
        X[prm.dim * i] = r * cos(aux);
        X[prm.dim * i + 1] = r * sin(aux);

        if (min_dist < R_exit / 100) {  // there is a body very close to i, so we remove it as well
          // set again randomly in the circle of radius R the body closest to i
          aux = 2 * M_PI * (double)rand() / RAND_MAX;
          r = sqrt(R * abs((double)rand() / RAND_MAX));
          X[prm.dim * closest] = r * cos(aux);
          X[prm.dim * closest + 1] = r * sin(aux);
        }
      }
    }

    // print progress
    if (numSteps % printSteps == 0) {
      cout << numSteps << " " << dt << " " << t << endl;
    }

    numSteps++;
  }

  return 0;
}
