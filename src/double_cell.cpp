#include "double_cell.h"

double_cell::double_cell (){
  x = y = delta_x = delta_y = -1.;
  //delta_x = delta_y = -1.;
}

/*double_cell::double_cell (int x, int y, int level, double delta_x, double delta_y){
  this->x = x;
  this->y = y;
  this->level = level;
  this->delta_x = delta_x;
  this->delta_y = delta_y;
  }*/

double_cell::double_cell (double x, double y, double dx, double dy){
  this->x = x;
  this->y = y;
  this->delta_x = dx;
  this->delta_y = dy;
}

double double_cell::get_double_cell_x(){
  return x;
}

double double_cell::get_double_cell_y(){
  return y;
}

double double_cell::get_double_cell_delta_x() {
  return delta_x;
}

double double_cell::get_double_cell_delta_y() {
  return delta_y;
}

/*double double_cell::get_double_cell_delta_y() {
  return delta_y;
  }*/

/*double double_cell::get_double_cell_delta_x() {
  return delta_x;
  }*/

void double_cell::set_double_cell_pointer_to_list(list<double_cell *>::iterator p){
  this->pointer_to_list = p;

}

list<double_cell *>::iterator double_cell::get_double_cell_pointer_to_list(){
  return this->pointer_to_list;

}

double_cell ** double_cell::split (){
  /*double_cell * cie, *cid, *cse, *csd;
  int newlevel = this->level + 1;
  double_cell ** V = (double_cell **) malloc (sizeof (double_cell *) * 4);
  cie = new double_cell(2 * (this->x), 2 * (this->y), newlevel);
  cid = new double_cell(2 * (this->x) + 1, 2 * (this->y), newlevel);
  cse = new double_cell(2 * (this->x), 2 * (this->y) + 1, newlevel);
  csd = new double_cell(2 * (this->x) + 1, 2 * (this->y) + 1, newlevel);
  V[0] = cie;
  V[1] = cid;
  V[2] = cse;
  V[3] = csd;
  return V;*/
  return NULL;
}

void double_cell::print_double_cell () {
  printf ("(%.7f, %.7f):%.7f %.7f\n", x, y, delta_x, delta_y);
}

void double_cell::set_double_cell_pointer_to_vector(int index){
  pointer_to_vector = index;
}

int double_cell::get_double_cell_pointer_to_vector(){
  return pointer_to_vector;
}

void double_cell::set_center_variable(double (* f) (double x, double y)){
  this->variable = (*f)(x, y);
}
