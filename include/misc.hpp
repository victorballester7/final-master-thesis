#ifndef MISC_HPP
#define MISC_HPP

#include <string>

using namespace std;

void saveData(int n, double *X, string filename_output, void *param);

pair<double, double> standard_normal();

void resetBody(double *X, int i, double R, void *param);

void removeBody(double **X, bool **isInside, int j, int *n, void *param);

void addBody(double **X, bool **isInside, int *n, double R_in, double C,
             void *param);

// computes the square of the velocity of the fluid at the position (x,y) due
// to the point vortices
double uv2(int n, double *X, double x, double y, void *param);

void Energy(int n, double t, double *X, string filename_output, void *param);
// compute the energy of the system in
void EnergyProf(int n, double *X, int N_R, double *E, double R_max,
                string filename_output, void *param);

void EnergyFlux(double dt, int N_R, double *E0, double *E1, double R_max,
                string filename_output);

void NumVortices(int n, double *X, int N_R, double R_max,
                 string filename_output, void *param);

#endif
