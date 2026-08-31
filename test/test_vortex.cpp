#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"

#define PI 3.1415926535897
#define Re 50
#define Sc 1
#define Xsi 150
#define xc 0.5
#define yc 0.5

using namespace std;

double u (double x, double y, double t){
  double r, vtc;
  r = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
  vtc = Re*Sc*(1.0 - exp(-(r*r)/(4.0*Sc*t)));
  return ((y - yc)*vtc/r);
}

double v (double x, double y, double t){
  double r, vtc;
  r = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
  vtc = Re*Sc*(1.0 - exp(-(r*r)/(4.0*Sc*t)));
  return (-(x - xc)*vtc/r);
}

double phi(double x, double y, double t){
  double d, r;
  d = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
  r = 0.25;
  //return (-4*PI*PI*cos(2*PI*x*t)*sin(2*PI*y*t));
  return (tanh(Xsi*(r-d)));
}

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
  int nxb = 32;
  int nyb = 32;

  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0;
  xend = yend = 1.;
  double t0 = 0.0;
  D = new dominio (xbegin, ybegin, xend, yend);
  
  mesh * M;
  
  /******create a base mesh BASE x BASE *********/
  //you can find the value for BASE at mesh.h file 
  M = new mesh(D, number_of_levels, nxb, nyb);
  /*********************************************/
    
  //M->get_hash_table()->print_information();

  int ct;
  vector <double> vmax;
  double dx, dy, xd, yd, dt, phis, erro;
  list <cell *> * l;
  list <cell *>::iterator it;
  double tempo;
  vector <cell *> * V;

  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();

  dt = 0.005;
  tempo = 1.0;
  
  M->initialize_var(&u, &v, &phi, tempo, t0);
  vmax = M->max_propriedades();
  erro = 0.0;
  //Refinamento inicial -
  ct = 0;
  int ns = 0;
  for(int k = 0; k < number_of_levels; k++){ 
    for (int i = 0; i < number_of_levels -1; i++) {
      //printf("%d\n", i);
      l = M->get_list_cell_by_level(i);
      
      dx = fabs(xend - xbegin) / (nxb * pow(2, i));
      dy = fabs(yend - ybegin) / (nyb * pow(2, i));
      
      it = l->begin();
      
      while (it != l->end()){
	xd = xbegin + (((*it)->get_cell_x() + 0.5) * dx);
	yd = ybegin + (((*it)->get_cell_y() + 0.5) * dy);
	
	phis = (*it)->get_cell_phi();
	erro = max(erro, sqrt(pow(phis - df(xd, yd, tempo),2)));  
	//refinement
	if (phis < 0.99 - k*0.05 && phis > - 0.99 + k*0.05){
	  ns++;
	  it = M->split(*it);
	}
	else
	  it++;
	
      }
      
      //M->initialize_var(&u, &v, &df, tempo, t0);
    }
    M->initialize_var(&u, &v, &phi, 0, 0);
  }
  
  M->print_silo(0, &P);
  //M->print(1000);

  
  /*  for (ct = 0; ct < 200; ct++){
    erro = 0.0;
    //refinement & merge
    for (int i = 0; i < number_of_levels; i++) {
      l = M->get_list_cell_by_level(i);
      
      dx = fabs(xend - xbegin) / (nxb * pow(2, i));
      dy = fabs(yend - ybegin) / (nyb * pow(2, i));

      it = l->begin();
      
      while (it != l->end()){
	xd = xbegin + (((*it)->get_cell_x() + 0.5) * dx);
	yd = ybegin + (((*it)->get_cell_y() + 0.5) * dy);
	
	phis = (*it)->get_cell_phi();
	erro = max(erro, sqrt(pow(phis - df(xd, yd, tempo),2)));  
	//refinement
	if ((phis <= 28. && phis >= 15.) || (phis >= -28. && phis <= -15.)){
	//if ((fabs(velu)/vmax[0] > 0.5 + i*0.1 && fabs(velu)/vmax[0] < 0.5 + (i+1)*0.1)  || (fabs(velv)/vmax[1] > 0.5 + i*0.1 && fabs(velv)/vmax[1] < 0.5 + (i+1)*0.1)){
	if (fabs(velu)/vmax[0] > 0.5 + i*0.1 || fabs(velv)/vmax[1] > 0.5 + i*0.1){

	  if (i < number_of_levels - 1 //se i não é o último nível, então pode refinar ){
	    it = M->split(*it);
	  }
	  else {
	    it++;
	  }
	  
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
    M->initialize_var(&u, &v, &df, tempo, t0);
    M->print(ct);
    
    printf("%f %f\n", tempo, erro);
    tempo += dt; 
    }*/
  
  M->get_hash_table()->print_information();
  cout << ns << endl;
  return 0;
}
