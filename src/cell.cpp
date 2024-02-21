#include "cell.h"

cell::cell (){
  x = y = level = index = -1, bc = 0;
  phi0 = phi = velv = velu = 0.0;
  cp = 0;
}


cell::cell (int x, int y, int level){
  this->x = x;
  this->y = y;
  this->level = level;
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

int cell::get_cell_bc(){
  return bc;
}

void cell::set_cell_index(int ivalue){
  index = ivalue;
}

void cell::set_cell_bc(int tbc){
  bc = tbc;
}

double cell::get_cell_velu(){
  return velu;
}

double cell::get_cell_velv(){
  return velv;
}

double cell::get_cell_phi(){
  return phi;
}

double cell::get_cell_phi0(){
  return phi0;
}

int cell::get_cell_with_particle(){
  return cp;
}

list <cell *> * cell::get_cell_nb(){
  return nb;
}

void cell::set_cell_nb(list <cell *> * vnb){
  nb = vnb;
}

list <weight> * cell::get_cell_wpoisson(){
  return wpoisson;
}

void cell::set_cell_wpoisson(list <weight> * lw){
  wpoisson = lw;
}

void cell::set_cell_velu(double uvalue){
  velu = uvalue;
}

void cell::set_cell_velv(double vvalue){
  velv = vvalue;
}

void cell::set_cell_phi(double phivalue){
  phi = phivalue;
}

void cell::set_cell_phi0(double phi0value){
  phi0 = phi0value;
}

void cell::set_cell_with_particle(int wparticle){
  cp = wparticle;
}

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

  cie = new cell(2 * (this->x), 2 * (this->y), newlevel);
  cid = new cell(2 * (this->x) + 1, 2 * (this->y), newlevel);
  cse = new cell(2 * (this->x), 2 * (this->y) + 1, newlevel);
  csd = new cell(2 * (this->x) + 1, 2 * (this->y) + 1, newlevel);

  V[0] = cie;
  V[1] = cid;
  V[2] = cse;
  V[3] = csd;
  return V;
}

void cell::print_cell () {
  printf ("%d %d %d %d\n", x, y, level, index);
}


