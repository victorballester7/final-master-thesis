#ifndef MISC_HPP
#define MISC_HPP

#include <string>

using namespace std;

void readData(int *n, double **X, int status, void *param);

void saveData(int n, double *X, string filename_output, string filename_status,
              double t, void *param, uint plot_count);

void EnergyTimes(int count_energy, double t, string filename_output);

pair<double, double> standard_normal();

void resetBody(double *X, int i, double R, void *param);

void removeBody(double **X, int j, int *n, void *param);

void addBody(double **X, bool **isInside, int *n, double R_in, double C,
             void *param);

void addVortices(double **X, int *n, int n_extra, double R_in, double eps,
                 void *param);
// computes the square of the velocity of the fluid at the position (x,y) due
// to the point vortices
double uv2(int n, double *X, double x, double y, void *param);

void Energy(int n, double t, double *X, string filename_output, void *param,
            uint plot_count);
// compute the energy of the system in
void EnergyProf(int n, double *X, int N_R, double *E, double dr,
                string filename_output, void *param, uint plot_count,
                bool save_data);

void EnergyFlux(double dt, int N_R, double *E0, double *E1, double dr,
                string filename_output, uint plot_count);

void NumVortices(int n, double *X, int N_R, double dr, string filename_output,
                 void *param, uint plot_count);

#endif
