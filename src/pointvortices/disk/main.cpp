#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "../../../include/misc.hpp"
#include "../../../include/pointvortices_field.h"
#include "../../../include/rk78.h"

#define TOL 1e-6                 // tolerance for the RK78 method
#define COLLISION (0.05 * R_out) // tolerance for the collision detection
#define MAX_VORTICES 2000        // maximum number of vortices
#define RK78 1

#define SIGN(x) ((x) > 0 ? 1 : -1)
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

using namespace std;

int main(void) {
  const string filename_input =
      "config/pointvortices/disk/input.txt"; // name of the input file to read
                                             // the parameters
  const string filename_output =
      "data/pointvortices/disk/positions/positions"; // name of the output file
                                                     // to write the positions
                                                     // of the point vortices
  const string filename_output_E =
      "data/pointvortices/disk/EnergyProf/Energy"; // name of the output file to
                                                   // write the positions of the
                                                   // point vortices
  const string filename_output_Eflux =
      "data/pointvortices/disk/EnergyFlux/Energy"; // name of the output file to
                                                   // write the positions of the
                                                   // point vortices
  const string filename_output_numvortices =
      "data/pointvortices/disk/NumVortices/NumVortices"; // name of the output
                                                         // file to write the
                                                         // number of point
                                                         // vortices
  const string filename_output_misc = "data/pointvortices/disk/misc.txt";
  const string filename_output_energyBal =
      "data/pointvortices/disk/energy_bal.txt";
  const string space = "    "; // space to print the parameters
  pointvortices_params prm;    // parameters of the system
  prm.dim = 3; // dimension of the space (2 space dimensions + 1 circulation
               // dimension)
  int n;       // number of bodies

  double R_in;      // radius of the circle in which the bodies are placed
  double R_in_aux;  // auxiliary radius of the circle in which the bodies are
                    // reinserted
  double R_out;     // radius of the circle in which the bodies are simulated
  double dt;        // time step
  uint totalSteps;  // total number of steps
  uint outSteps;    // number of steps to write the positions of the bodies
  uint energySteps; // number of steps to write the energy profile
  uint printSteps;  // number of steps to print the progress
  uint inputSteps = 50; // frequency of the input steps

  int N_R = 50;  // number of points to discretize the radial direction
  double *E;     // energy profile
  double *E_aux; // energy profile at time t
  bool *isInner; // boolean array to check if the body is inside the inner
                 // cercle circle

  E = (double *)malloc(N_R * sizeof(double));
  E_aux = (double *)malloc(N_R * sizeof(double));

  time_t t0 = time(0);

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
    file_input >> tmp >> prm.alpha;
    file_input >> tmp >> R_in;
    file_input >> tmp >> R_out;
    file_input >> tmp >> dt;
  }
  file_input.close();

  file_output.open(filename_output_misc);
  if (file_output.is_open()) {
    file_output << R_out << " " << outSteps << endl;
  }
  file_output.close();

  // print parameters
  // ------------- Print plot setup -----------------
  cout << "Number of steps:    " << space << totalSteps << endl;
  cout << "Number of votices:  " << space << n << endl;
  cout << "Initial radius:     " << space << R_in << endl;
  cout << "Exit radius:        " << space << R_out << endl;
  cout << "Drag coefficient:   " << space << prm.alpha << endl;

  isInner = (bool *)malloc(n * sizeof(bool));

  int N = prm.dim * n; // dimension of the field
  prm.EPS = 1e-3;      // softening parameter
  // -------------------------------------------

  double t = 0.0;
  uint numSteps = 0;

  // allocate memory
  double *X; // vector of positions (x1,y1,C1,x2,y2,C2,...,xn,yn,Cn)
  X = (double *)malloc(N * sizeof(double));

  double dist;

  // set circulations
  double eps = 1.; // small perturbation to avoid the same circulations
  pair<double, double> z;
  // seed the random number generator
  srand(time(NULL));
  double aux, r;

  // set initial conditions (random in the circle of radius R)
  for (int i = 0; i < n; i++) {
    if (i % 2 == 0) {
      z = standard_normal();
      X[prm.dim * i + 2] = eps * z.first;
      X[prm.dim * (i + 1) + 2] = -eps * z.first;
    }
    aux = 2 * M_PI * (double)rand() / RAND_MAX;
    r = R_in * sqrt(abs((double)rand() /
                        RAND_MAX)); // random radius. We do not do it uniformly
                                    // to avoid having more points per area in
                                    // the smaller subradii. With the sqrt we
                                    // have more points in the outer layers.
    X[prm.dim * i] = r * cos(aux);  // x position
    X[prm.dim * i + 1] = r * sin(aux); // y position
    isInner[i] = true;
  }

  double hmin = dt / 100;
  double hmax = dt * 100;

  int n0 = n;
  int n_extra = 2; // should be even

  while (numSteps < totalSteps) {
    // save data
    if (numSteps % outSteps == 0) {
      saveData(n, X, filename_output, &prm);
    }

    // plot energy profile
    if (numSteps % energySteps == 0) {
      EnergyProf(n, X, N_R, E, R_out, filename_output_E, &prm);
      Energy(n, t, X, filename_output_energyBal, &prm);
      NumVortices(n, X, N_R, R_out, filename_output_numvortices, &prm);
    }

    if ((numSteps - 1) % energySteps == 0) {
      EnergyProf(n, X, N_R, E_aux, R_out, filename_output_E, &prm);
      EnergyFlux(dt, N_R, E, E_aux, R_out, filename_output_Eflux);
    }

    // simulate
#ifdef RK78
    if (rk78(&t, X, &dt, hmin, hmax, TOL, N, pointvortices_field, &prm) != 0) {
      cout << "Error in the integration" << endl;
      free(X);
      free(E);
      free(E_aux);
      free(isInner);
      return 1;
    }
#else
    if (rk4(&t, X, &dt, N, pointvortices_field, &prm) != 0) {
      cout << "Error in the integration" << endl;
      free(X);
      free(E);
      free(E_aux);
      free(isInner);
      return 1;
    }
#endif

    // check if a body has exited the domain
    if (numSteps % energySteps != 0) {
      for (int i = 0; i < n; i++) {
        if (X[prm.dim * i] * X[prm.dim * i] +
                X[prm.dim * i + 1] * X[prm.dim * i + 1] >
            R_out * R_out) {
          // search the closest point to i-th body (advection of 2 bodies in a
          // straight line)

          // min_dist = R_exit; // Radius of the circle
          // closest = -1;
          for (int j = i + 1; j < n; j++) {
            dist = (X[prm.dim * i] - X[prm.dim * j]) *
                       (X[prm.dim * i] - X[prm.dim * j]) +
                   (X[prm.dim * i + 1] - X[prm.dim * j + 1]) *
                       (X[prm.dim * i + 1] - X[prm.dim * j + 1]);

            if (dist < COLLISION) {
              // cout << "adeuuuuuuuuuuuuuu1" << endl;
              removeBody(&X, &isInner, j, &n, &prm);
              // cout << "adeuuuuuuuuuuuuuu2" << endl;
              j--;
            }
          }
          // set again randomly in the circle of radius R the body i
          // cout << "mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmprova1" <<
          // endl;
          removeBody(&X, &isInner, i, &n, &prm);
          // cout << "mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmprova2" <<
          // endl;
          i--;
        }
      }
    }
    if (numSteps % inputSteps == 0 && numSteps > 0) {
      // cout << "hello1" << endl;
      addVortices(&X, &n, n_extra, R_in, eps, &prm);
      // cout << "hello2" << endl;
    }
    // print progress
    if (numSteps % printSteps == 0) {
      cout << numSteps << " " << n << " " << n0 << " " << dt << " " << t
           << endl;
    }

    numSteps++;
    if (n > MAX_VORTICES) {
      cout << "Maximum number of vortices reached (" << n << ")" << endl;
      break;
    }
    N = prm.dim * n;
  }

  time_t t1 = time(0);
  cout << "Time employed: " << difftime(t1, t0) << " seconds" << endl;

  // free memory
  free(X);
  free(E);
  free(E_aux);
  free(isInner);

  return 0;
}
