#include "cell.h"

cell::cell (){
  x = y = level = -1;
  //delta_x = delta_y = -1.;
}

/*cell::cell (int x, int y, int level, double delta_x, double delta_y){
  this->x = x;
  this->y = y;
  this->level = level;
  this->delta_x = delta_x;
  this->delta_y = delta_y;
  }*/

cell::cell (int x, int y, char kind[], int level){
  this->x = x;
  this->y = y;
  this->kind[0] = kind[0];
  this->kind[1] = kind[1];
  this->last_kind[0] = kind[0];
  this->last_kind[1] = kind[1];
  this->level = level;
}

/*cell::cell (int x, int y, char kind[], char last_kind[], int level){
  this->x = x;
  this->y = y;
  this->kind[0] = kind[0];
  this->kind[1] = kind[1];
  this->last_kind[0] = last_kind[0];
  this->last_kind[1] = last_kind[1];
  this->level = level;
  }*/

int cell::get_cell_x(){
  return x;
}

void cell::set_cell_kind(char kind []){
  this->kind[0] = kind[0];
  this->kind[1] = kind[1];
}

int cell::get_cell_y(){
  return y;
}

char * cell::get_cell_kind(){
  return kind;
}

char * cell::get_cell_last_kind(){
  return last_kind;
}

int cell::get_cell_level() {
  return level;
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
  cell ** V = (cell **) malloc (sizeof (cell *) * 4);
  char kind [2];
  
  kind[0] = 'b'; kind[1] = 'l';
  cie = new cell(2 * (this->x), 2 * (this->y), kind, newlevel);
  kind[0] = 'b'; kind[1] = 'r';
  cid = new cell(2 * (this->x) + 1, 2 * (this->y), kind, newlevel);
  kind[0] = 't'; kind[1] = 'l';
  cse = new cell(2 * (this->x), 2 * (this->y) + 1, kind, newlevel);
  kind[0] = 't'; kind[1] = 'r';
  csd = new cell(2 * (this->x) + 1, 2 * (this->y) + 1, kind, newlevel);
  V[0] = cie;
  V[1] = cid;
  V[2] = cse;
  V[3] = csd;
  return V;
}

void cell::print_cell () {
  printf ("(%d, %d):%d - %c%c\n", x, y, level, kind[0], kind[1]);
}


