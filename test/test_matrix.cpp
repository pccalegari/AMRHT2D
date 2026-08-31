#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"
#include "mls.h"
#include "ode.h"

#define PI 3.1415926535897

using namespace std;

double phi(double x, double y, double t){
 
  //return (exp(-t)*cos(2*PI*x)*sin(2*PI*y));
  //return (0.5*(1.0 + tanh(Xsi*(y - yc))));
  return cos(2*PI*x)*sin(2*PI*y);
}

double f (double x, double y) { 
  return 8.0*cos(2 * PI * x) * sin(2 * PI * y);
}

int main (){
  
  int number_of_levels = 4;
  int nxb = 8;
  int nyb = 8;

  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = 0.0;
  xend = 1.0;
  ybegin = 0.0;
  yend = 1.0;
  D = new dominio (xbegin, ybegin, xend, yend);
  
  mesh * M;
  
  /******create a base mesh BASE x BASE *********/
  //you can find the value for BASE at mesh.h file 
  M = new mesh(D, number_of_levels, nxb, nyb);
  /*********************************************/
  
  vector <double> vmax, xv, yv, uv, vv;
  list <cell *> * l;
  list <cell *>::iterator it, itv, itb;
  vector <cell *> *V;
    
  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();
  double tempo = 1.0, yd, xd, dx, dy;
  int ix, iy, ncell, ncellv;

  for (int i = 0; i < number_of_levels - 1; i++) {
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));

    it = l->begin();
    
    while (it != l->end()){
      ix = (*it)->get_cell_x();
      iy = (*it)->get_cell_y();
      xd = xbegin + ix * dx;
      yd = ybegin + iy * dy;
      
      //refinement
      if ((xd >= 0.75 + dx*i && yd < 0.5 - dy*i) || (xd >= 0.5 + dx*i && xd < 0.75 - dx*i && yd >= 0.25 + dy*i && yd < 0.5 - dy*i) ){
	
	it = M->split(*it);
      }
      else {
	it++;
      }
      //***
    }
  }

  ncell = M->counting_mesh_cell();
  cout << ncell << endl;
  cout << nxb << endl;
  M->print_mesh();
  ncellv = M->neighbours_all_cell();

  M->create_unstructured_mesh(&phi, tempo);

  return 0;
}


