#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "ode.h"

using namespace std;

ode::ode(){
  tini = sini = 0.0;
  tend = 0.0;
  //double (*sfn)(double t, double phi);  // Funcao fonte
}

ode::ode(double t0, double phi0, double tn, double (*f)(double, double)){ // Construtor
  this->tini = t0;
  this->sini = phi0;
  this->tend = tn;
  this->sfn = f;
}

double ode::rk1 (double dt){
  double phi, phi0, t;
  phi = -1.0;
  phi0 = sini;
  t  = tini;
  while(t != tend){
    phi = phi0 + dt*sfn(t, phi0);
    phi0 = phi;
    t += dt;
    dt = min(dt, tend - t);
  }
  //printf("%f %f\n", t, phi);
  return phi;
}

double ode::rk2M (double dt){
  double phi, phi0, t, k1;
  phi = -1.0;
  phi0 = sini;
  t  = tini;
  while(t != tend){
    k1 = sfn(t, phi0);
    phi = phi0 + dt*sfn(t + 0.5*dt, phi0 + 0.5*dt*k1);
    phi0 = phi;
    t += dt;
    dt = min(dt, tend - t);
  }
  //printf("%f %f\n", t, phi);
  return phi;
}

double ode::rk2A (double dt){
  double k1, k2, phi, phi0, t;
  phi = -1.0;
  phi0 = sini;
  t  = tini;
  while(t != tend){
    k1 = sfn(t, phi0);
    k2 = sfn(t + dt, phi0 + dt*k1);
    phi = phi0 + 0.5*dt*(k1 + k2);
    phi0 = phi;
    t += dt;
    dt = min(dt, tend - t);
  }
  //printf("%f %f\n", t, phi);
  return phi;
}

double ode::rk4 (double dt){
  double k1, k2, k3, k4, phi, phi0, t;
  phi = -1.0;
  phi0 = sini;
  t = tini;
  while(t != tend){
    k1 = dt*sfn(t, phi0);
    k2 = dt*sfn(t + 0.5*dt, phi0 + 0.5*k1);
    k3 = dt*sfn(t + 0.5*dt, phi0 + 0.5*k2);
    k4 = dt*sfn(t + dt, phi0 + k3);
    phi = phi0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    phi0 = phi;
    //printf("%f %f %f\n", dt, t, phi);
    t += dt;
    dt = min(dt, tend - t);
  }
  //printf("%f %f\n", t, phi);
  return phi;
}
  
