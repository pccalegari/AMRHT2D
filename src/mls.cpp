#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mls.h"
#include "point.h"

using namespace std;

double delta(double r, int op){
  double deltas = 0.0;
  if(op == 1)
    deltas = delta_s1(r);
  else if(op == 2)
    deltas = delta_s2(r);
  else if(op == 3)
    deltas = delta_s3(r);
  return deltas;
}

double delta_s2(double r){
  double delta = 0.0;
  r = fabs(r);
  if(r >= 0.0 && r < 1.0)
    delta = (1.0/8.0)*(3.0 - 2.0*r + sqrt(1.0 + 4.0*r - 4.0*r*r));
  else if(r >= 1.0 && r <= 2.0)
    delta = (1.0/8.0)*(5.0 - 2.0*r - sqrt(12.0*r - 7.0 - 4.0*r*r));
  return delta;
}

double delta_s1(double r){
  double delta = 0.0;
  r = fabs(r);
  if(r >= 0 && r <= 0.5)
    delta = 2/3 - 4*r*r + 4*r*r*r;
  else if(r > 0.5 && r <= 1.0)
    delta = 4/3 - 4*r + 4*r*r - (4/3)*r*r*r;
  return delta;
}

double delta_s3(double r){
  double delta = 0.0;
  r = fabs(r);
  if((r >= 0) && (r <= 1))
    delta = 1 - r;
  return delta;
}

double delta_inter(double x, double y, vector <double> xv, vector <double> yv, vector <double> phi, double d, int op){
  int ndim = phi.size();
  double aphi = 0.0, deltax, deltay;
  for(int i = 0; i < ndim; i++){
    deltax = delta((xv[i] - x)/d,op);
    deltay = delta((yv[i] - y)/d,op);
    aphi += phi[i]*deltax*deltay;
  }
  
  return aphi;
}

double max_var(vector <double> var){
  int ndim = var.size();
  double maxvar = 0.0;
  for(int i = 0; i < ndim; i++){
    if(maxvar < fabs(var[i]))
       maxvar = fabs(var[i]);
  }
  
  return maxvar;
}

/*void mls::pesos_inter_diag(point p, cell c){
  
  q = new point(0,0);
  double h1 = dist(p, q);
  //double h2 = dist(p, r);
}

*/  
  
