#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

double rk1(double dt, double phi0, double k1);
double rk2(double dt, double phi0, double k1, double k2);
double rk4(double dt, double phi0, double k1, double k2, double k3, double k4);

