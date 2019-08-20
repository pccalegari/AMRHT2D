#include "cell.h"

cell::cell (){
  x = y = level = index = -1;
  //delta_x = delta_y = -1.;
}

/*cell::cell (int x, int y, int level, double delta_x, double delta_y){
  this->x = x;
  this->y = y;
  this->level = level;
  this->delta_x = delta_x;
  this->delta_y = delta_y;
  }*/

cell::cell (int x, int y, int level, int index){
  this->x = x;
  this->y = y;
  this->level = level;
  this->index = index;
}

int cell::get_cell_x(){
  return x;
}

int cell::get_cell_y(){
  return y;
}

int cell::get_cell_level() {
  return level;
}

int cell::get_cell_index(){
  return index;
}

void cell::put_cell_index(int ivalue){
  index = ivalue;
}
  
/*double cell::get_cell_delta_y() {
  return delta_y;
  }*/

/*double cell::get_cell_delta_x() {
  return delta_x;
  }*/

void cell::set_cell_pointer_to_list(list<cell *>::iterator p){
  pointer_to_list = p;

}

list<cell *>::iterator cell::get_cell_pointer_to_list(){
  return pointer_to_list;

}

cell ** cell::split (){
  cell * cie, *cid, *cse, *csd;
  int newlevel = this->level + 1;
  int index_p = this->index;
  //double new_delta_x = this->delta_x/2.;
  //double new_delta_y = this->delta_y/2.;
  cell ** V = (cell **) malloc (sizeof (cell *) * 4);
  cie = new cell(2 * (this->x), 2 * (this->y), newlevel, index_p);
  cid = new cell(2 * (this->x) + 1, 2 * (this->y), newlevel, index_p);
  cse = new cell(2 * (this->x), 2 * (this->y) + 1, newlevel, index_p);
  csd = new cell(2 * (this->x) + 1, 2 * (this->y) + 1, newlevel, index_p);
  V[0] = cie;
  V[1] = cid;
  V[2] = cse;
  V[3] = csd;
  return V;
}

void cell::print_cell () {
  printf ("(%d, %d):%d %d\n", x, y, level, index);
}


