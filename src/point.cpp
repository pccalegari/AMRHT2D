#include "point.h"
#include <cmath>

point::point (){
  x = y = 0.0;
}


point::point (double x, double y){
  this->x = x;
  this->y = y;
}

double point::get_point_x(){
  return x;
}

double point::get_point_y(){
  return y;
}

void point::set_point_x(double xp){
  x = xp;
}

void point::set_point_y(double yp){
  y = yp;
}

double point::dist(point p, point q){
  double px, py, qx, qy;
  px = p.get_point_x();
  py = p.get_point_y();
  qx = q.get_point_x();
  qy = q.get_point_y();
  double d = sqrt((px - qx)*(px - qx) + (py - qy)*(py - qy));
  return d;
}


