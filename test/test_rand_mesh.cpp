#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"

#define PI 3.1415926535897

using namespace std;

double f (double x, double y) { 
  return cos(2 * PI * x) * sin(2 * PI * y);
}

double df (double x, double y, double tempo) {
  double xc, yc, r, phi, d;
  xc = 0.5;
  yc = 0.5;
  r = 0.25;
  d = sqrt((x-xc)*(x-xc) +(y-yc)*(y-yc));
  phi = tanh(75*(r-d));
  //return -4 * PI * PI * f(tempo * x, tempo * y);
    return(phi);
}

int main (){

  int number_of_levels = 4;
  int nxb = 16;
  int nyb = 16;
  cell * c;
  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.;
  xend = yend = 1.;
  D = new dominio (xbegin, ybegin, xend, yend);
  printf("%d %d %d\n", nxb, nyb, number_of_levels);
  mesh * M;
  
  /******create a base mesh BASE x BASE *********/
  //you can find the value for BASE at mesh.h file 
  M = new mesh(D, number_of_levels, nxb, nyb);
  /*********************************************/
    
  //M->get_hash_table()->print_information();

  
  double dx, dy, x, y;
  list <cell *> * l;
  list <cell *>::iterator it;
  
  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();

  //srand(time(NULL));
  srand(12345);
  int r;
  
  double dxf = fabs(xend - xbegin) / (nxb * pow(2, number_of_levels - 1));
  for (int i = 0; i < number_of_levels -1; i++) {
    
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    
    it = l->begin();
    
    while (it != l->end()){
      r = rand();
      //xd = xbegin + (((*it)->get_cell_x()) * dx);
      //yd = ybegin + (((*it)->get_cell_y()) * dy);

      if (r % 3 == 0 && r % 2 == 1){
	c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
	assert(c != NULL);
	it = M->split(c);
      }
      else
	it++;
      
    }
    //M->initialize_var(&u, &v, &df, tempo, t0);
  }
  
  M->create_unstructured_mesh(&df, dxf);
  //M->get_hash_table()->print_information();
  int ncell = M->counting_mesh_cell();
  printf("%d\n", ncell);

  for (int i = 0; i < number_of_levels; i++) {
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    
    for (list <cell *>::reverse_iterator it = l->rbegin(); it != l->rend(); it++) {
      x = xbegin + (*it)->get_cell_x() * dx;
      y = ybegin + (*it)->get_cell_y() * dy;
      printf("%f %f %d\n",x, y, (*it)->get_cell_level());
      
    }
  }

  
  return 0;
}
