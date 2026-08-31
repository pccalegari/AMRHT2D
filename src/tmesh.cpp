#include <cmath>
#include <cassert>
#include "mesh.h"
#include <vector>
#include "cell.h"

tmesh::vertice(){
  vx = vy = 0.0;
  cell = NULL;
  index = -1;
}

tmesh::aresta(){
  vert1 = NULL;
  vert2 = NULL;
  index = -1;
}

tmesh::elemento(){
  vert1 = NULL;
  vert2 = NULL;
  vert3 = NULL;
  index = -1;
}

tmesh::vertice(double x, double y, int id, cell * c){
  this->vx = x;
  this->vy = y;
  this->cell = c;
  this->index = id;
}

double tmesh::get_vx(){
  return x;
}

double tmesh::get_vy(){
  return y;
}

int tmesh::get_vid(){
  return id;
}

cell * tmesh::get_vcell() {
  return c;
}

tmesh::aresta(){
}
