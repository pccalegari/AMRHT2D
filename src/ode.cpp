#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "ode.h"

using namespace std;

// rk1 ou rk2M k1 = f(t,phi) k2 = f(t+0.5*dt,phi+0.5*dt*k1) 
double rk1 (double dt, double phi0, double k1){
  double phi;
  phi = phi0 + dt*k1;
  return phi;
}

double rk2 (double dt, double phi0, double k1, double k2){
  double phi;
  phi = phi0 + 0.5*dt*(k1 + k2);
  return phi;
}

double rk4 (double dt, double phi0, double k1, double k2, double k3, double k4){
  double phi;
  phi = phi0 + dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
  return phi;
}
  
