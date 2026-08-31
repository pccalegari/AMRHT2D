#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"
#include "particle.h"
#include "print.h"

#define PI 3.1415926535897

using namespace std;

double f (double x, double y) { 
  return cos(2 * PI * x) * sin(2 * PI * y);
}

double df (double x, double y, double tempo) {
  return -4 * PI * PI * f(tempo * x, tempo * y);
}

int main (){
  
  list <particle *> P;

  particle * p = new particle();

  P.push_back(p);

  particle * q = P.back();

  q->print_particle();
  
  int number_of_levels = 4;
  int nxb = 64;
  int nyb = 64;

  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = -1.;
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
  
  for (tempo = -1; tempo <= 0; tempo += 0.005){
    //refinement & merge
    for (int i = number_of_levels - 1; i >= 0 ; i--) {
      l = M->get_list_cell_by_level(i);
      
      dx = fabs(xend - xbegin) / (nxb * pow(2, i));
      dy = fabs(yend - ybegin) / (nyb * pow(2, i));

      it = l->begin();
      
      while (it != l->end()){
	
	xd = xbegin + ((*it)->get_cell_x() * dx);
	yd = ybegin + ((*it)->get_cell_y() * dy);
	
	//refinement
	if ((df(xd + (dx / 2.), yd + (dy / 2.), tempo) <= 25 && df(xd + (dx / 2.), yd + (dy / 2.), tempo) >= 18) ||
	    (df(xd + (dx / 2.), yd + (dy / 2.), tempo) >= -25 && df(xd + (dx / 2.), yd + (dy / 2.), tempo) <= -18)){
	  if (i < number_of_levels - 1 /*se i não é o último nível, então pode refinar*/){
	    it = M->split(*it);
	  }
	  else {
	    it++;
	  }
	  //***
	}
	//merge
	else{
	  if (i > 0) {//se i não é o primeiro nível, então pode fazer merge
	    V = M->siblings((*it));
	    if (V != NULL){
	      it = M->merge(V, *it);
	      V->clear();
	    }
	    else it++;
	  }
	  else it++;
	}
      }
    }
    M->create_unstructured_mesh(&df, tempo);
    //M->get_hash_table()->print_information();
  }
  
  M->create_unstructured_mesh(&df, tempo);
  //M->get_hash_table()->print_information();
  
  return 0;
}
