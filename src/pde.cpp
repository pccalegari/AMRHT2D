#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "dominio.h"
#include "mesh.h"

using namespace std;

pde::pde(){
}

pde::pde(mesh * M, double (*fx1)(double, double, double), double (*fx2)(double, double, double), double (*fy1)(double, double, double), double (*fy2)(double, double, double)){ // Construtor
  this->M = M;
  this->bc_south = fx1;
  this->bc_north = fx2;
  this->bc_west = fy1;
  this->bc_east = fy2;
}

double pde::dbc_phi(cell * c, double t){
  double x, y, xbegin, xend, ybegin, yend, dx, dy;
  int l, i, j;
  double value = 1.0e+10, phiv;
  xbegin = get_dominio()->get_xbegin();
  xend = get_dominio()->get_xend();
  ybegin = get_dominio()->get_ybegin();
  yend = get_dominio()->get_yend();
  l = (*c)->get_cell_level();
  dx = fabs(xend - xbegin) / (nxb * pow(2, l));
  dy = fabs(yend - ybegin) / (nyb * pow(2, l));
  i = (*c)->get_cell_x();
  j = (*c)->get_cell_y();
  phiv = (*c)->get_cell_phi();
  
  if(i == 0){
    x = xbegin;
    y = ybegin + 0.5 * j * dy;
    value = 2.0*bc_west(x, y, t) - phiv;
  }
  else if(i == pow(2,l)*nxb - 1){
    x = xend;
    y = ybegin + 0.5 * j * dy;
    value = 2.0*bc_east(x, y, t) - phiv;
  }
  else if(j == 0){
    x = xbegin + 0.5 * i * dx;
    y = ybegin;
    value = 2.0*bc_south(x, y, t) - phiv;
  }
  else if(j == pow(2,l)*nxb - 1){
    x = xbegin + 0.5 * i * dx;
    y = yend;
    value = 2.0*bc_north(x, y, t) - phiv;
  }
  return value
}

double pde::nbc_phi(cell * c, double t){
  double x, y, xbegin, xend, ybegin, yend, dx, dy;
  int l, i, j;
  double value = 1.0e+10, phiv;
  xbegin = get_dominio()->get_xbegin();
  xend = get_dominio()->get_xend();
  ybegin = get_dominio()->get_ybegin();
  yend = get_dominio()->get_yend();
  l = (*c)->get_cell_level();
  dx = fabs(xend - xbegin) / (nxb * pow(2, l));
  dy = fabs(yend - ybegin) / (nyb * pow(2, l));
  i = (*c)->get_cell_x();
  j = (*c)->get_cell_y();
  phiv = (*c)->get_cell_phi();
  
  if(i == 0){
    x = xbegin;
    y = ybegin + 0.5 * j * dy;
    value = phiv - dx*bc_west(x, y, t);
  }
  else if(i == pow(2,l)*nxb - 1){
    x = xend;
    y = ybegin + 0.5 * j * dy;
    value = phiv + dx*bc_east(x, y, t);
  }
  else if(j == 0){
    x = xbegin + 0.5 * i * dx;
    y = ybegin;
    value = phiv - dy*bc_south(x, y, t);    
  }
  else if(j == pow(2,l)*nxb - 1){
    x = xbegin + 0.5 * i * dx;
    y = yend;
    value = phiv + dy*bc_north(x, y, t);
  }
  return value
}

  
