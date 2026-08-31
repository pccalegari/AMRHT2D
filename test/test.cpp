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

  int number_of_levels = 3;
  int nxb = 32;
  int nyb = 32;

  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.;
  xend = yend = 1.;
  D = new dominio (xbegin, ybegin, xend, yend);
  
  mesh * M;
  
  /******create a base mesh BASE x BASE *********/
  //you can find the value for BASE at mesh.h file 
  M = new mesh(D, number_of_levels, nxb, nyb);
  /*********************************************/
    
  //M->get_hash_table()->print_information();

  
  double dx, dy, xd, yd;
  list <cell *> * l;
  list <cell *>::iterator it;
  double tempo;
  vector <cell *> * V;

  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();

  double dxf = fabs(xend - xbegin) / (nxb * pow(2, number_of_levels - 1));
  for (int i = 0; i < number_of_levels -1; i++) {
    
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    
    it = l->begin();
    
    while (it != l->end()){
      xd = xbegin + (((*it)->get_cell_x()) * dx);
      yd = ybegin + (((*it)->get_cell_y()) * dy);

      if (df(xd+0.5*dx, yd+0.5*dy, 0.0) < 0.999 - i*0.05 && df(xd+0.5*dx, yd+0.5*dy, 0.0) > - 0.999 + i*0.05){
      //if (df(xd, yd, 0.0) < 0.9 && df(xd, yd, 0.0) > - 0.9){
      	it = M->split(*it);
      }
      else
	it++;
      
    }
    //M->initialize_var(&u, &v, &df, tempo, t0);
  }
  
  M->create_unstructured_mesh(&df, dxf);
  M->get_hash_table()->print_information();
  
  return 0;
}
