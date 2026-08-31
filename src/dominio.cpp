#include "dominio.h"

dominio::dominio (){
  xbegin = ybegin = xend = yend = 0.;
}

dominio::dominio (double xb, double yb, double xe, double ye){
  xbegin = xb;
  ybegin = yb;
  xend = xe;
  yend = ye;
}

void dominio::set_xbegin(double x){
  xbegin = x;
}

void dominio::set_ybegin(double y){
  ybegin = y;
}

void dominio::set_xend(double x){
  xend = x;
}

void dominio::set_yend(double y){
  yend = y;
}

double dominio::get_xbegin(){
  return xbegin;
}

double dominio::get_ybegin(){
  return ybegin;
}

double dominio::get_xend(){
  return xend;
}

double dominio::get_yend(){
  return yend;
}
