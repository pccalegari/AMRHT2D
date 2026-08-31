#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"
#include "mls.h"
#include "ode.h"

#define PI 3.1415926535897
#define Re 50.0
#define Sc 10.
#define Xsi 200.0 //200.0
#define xc 0.
#define yc 0.
#define Ca 24.0
#define rdens 1.18/2470.0

using namespace std;


double u (double x, double y, double t){
  double r, vtc;
  r = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
  vtc = (Re*Sc/(2*PI*r))*(1.0 - exp(-(r*r)/(4.0*Sc*(t+1.0e-6))));
  return ((y - yc)*vtc/r);
  //return(-r*sin(t));
}

double v (double x, double y, double t){
  double r, vtc;
  r = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc));
  vtc = (Re*Sc/(2*PI*r))*(1.0 - exp(-(r*r)/(4.0*Sc*(t+1.e-6))));
  return (-(x - xc)*vtc/r);
  //return(r*cos(t));
}

double phi(double x, double y, double t){
 
  //return (exp(-t)*cos(2*PI*x)*sin(2*PI*y));
  //return (0.5*(1.0 + tanh(Xsi*(y - yc))));
  return 1.0;
}

double f (double x, double y) { 
  return cos(2 * PI * x) * sin(2 * PI * y);
}

double df (double x, double y, double t) {
  return (-exp(-t)*cos(2*PI*x)*sin(2*PI*y));
}

int main (){
  clock_t start0, finish0, start, finish, start1, finish1;

  start = clock();
  list <particle *> P;
  
  int npart = 1000;
  //particle * particle_bck;

  int number_of_levels = 4;
  int nxb = 16;
  int nyb = 16;

  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = -0.5;
  xend = 0.5;
  ybegin = -0.5;
  yend = 0.5;
  D = new dominio (xbegin, ybegin, xend, yend);
    
  mesh * M;
  
  /******create a base mesh BASE x BASE *********/
  //you can find the value for BASE at mesh.h file
  start0 = clock();
  M = new mesh(D, number_of_levels, nxb, nyb);
  finish0 = clock();
  
  /*********************************************/
  
  srand(123456);
  double raio;
  double xp, yp;//, up, vp, ufp, vfp, taup, nvelup;
  
  vector <double> vmax, xv, yv, uv, vv;
  double dt; //dx, dy, dt;//, velu, velv;
  double k1u, k1v, k2u, k2v, k3u, k3v, k4u, k4v;
  list <cell *> * l;
  list <cell *>::iterator it, itv, itb;
  double tempo = 0.0, tf = 1.0;
  //list <cell *> * lv, ln;
  //cell ** sv;
  vector <cell *> *V;
  int ctmax = 300;
  
  cell * cparticle, * cparticle_bck;
 
  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();

  raio = fabs((xend - xbegin))/2.0 - 0.05;
    
  for(int k = 0; k < npart; k++){
    particle * p = new particle();
    if(rand()%2 == 0){
      p->set_particle_x(raio*cos(rand())+xc);
      p->set_particle_y(yc+0.001*sin(rand()));
    }
    else{
      p->set_particle_x(raio*cos(rand())+xc);
      p->set_particle_y(yc-0.001*sin(rand()));
    }
    p->set_particle_radius(20.0e-6);
    xp = p->get_particle_x();
    yp = p->get_particle_y();
    p->set_particle_vx(u(xp,yp,tempo));
    p->set_particle_vy(v(xp,yp,tempo));
    
    //refina a celula que contem a particula ate o nivel mais fino
    for(int k = 0; k < number_of_levels - 1; k++){
      cparticle = M->search_particle_cell(p, k);
      if(cparticle != NULL){
	M->split(cparticle);
      }
    }
    //Neste ponto a particula esta no nivel mais fino
    p->set_particle_level(number_of_levels - 1);
    cparticle = M->search_particle_cell(p, p->get_particle_level());
    cparticle->set_cell_with_particle(1);
    /*p->print_particle();
    lv = M->neighbours(cparticle);
    for(itv = lv->begin(); itv != lv->end(); itv++){
      (*itv)->print_cell();
      if((*itv)->get_cell_level() != number_of_levels - 1){
	ln.push_front(*itv);
      }
    }
    while(ln.size() != 0){
      itv = ln.begin();
      ln.pop_front();
      sv = M->split_return_new_cells(*itv);
      
      for(int k = 0; k < 4; k++){
	if(sv[k]->get_cell_level() != number_of_levels - 1)
	  ln.push_front(sv[k]);
      }
      }*/
    
    P.push_back(p);
  }

  int ct = 0;
  double cputime = 0.0;
  M->initialize_var(&u, &v, &phi, tempo, tempo);
  //printf("%d %.8e\n", ct, tempo);
  //M->print_silo(ct, &P);
  
  double h = min(fabs(xend - xbegin) / (nxb * pow(2, number_of_levels-1)), fabs(yend - ybegin) / (nyb * pow(2, number_of_levels-1)));
  vmax = M->max_propriedades();
  dt = 0.5*h/max(vmax[0],vmax[1]);
  
  while(tempo < tf && ct <= ctmax){
    //printf("%d %E\n", ct+1, tempo);
  
    //Localiza particula e marca posicao da particulas e vizinhas
    // -> "deve ser uma funcao" 
    l = M->get_list_cell_by_level(number_of_levels - 1);
    it = l->begin();
    
    while (it != l->end()){
      (*it)->set_cell_with_particle(0);
      it++;
    }

    
    for(list <particle *>::iterator itp = P.begin(); itp != P.end(); itp++){
      
      xp = (*itp)->get_particle_x();
      yp = (*itp)->get_particle_y();
      /*cparticle_bck = M->search_particle_cell((*itp), number_of_levels - 1);
      cparticle_bck->set_cell_with_particle(1);
      particle_bck = new particle(xp, yp, 0, 0, (*itp)->get_particle_level());
      Pbck.push_back(particle_bck);*/
      k1u = u(xp,yp,tempo);
      k1v = v(xp,yp,tempo);
      // RK1
      //(*itp)->set_particle_x(rk1(dt, xp, k1u));
      //(*itp)->set_particle_y(rk1(dt, yp, k1v));
      // RK2
      //k2u = u(xp + dt*k1u, yp + dt*k1v, tempo + dt);
      //k2v = v(xp + dt*k1u, yp + dt*k1v, tempo + dt);
      //(*itp)->set_particle_x(rk2(dt, xp, k1u, k2u));
      //(*itp)->set_particle_y(rk2(dt, yp, k1v, k2v));
      k2u = u(xp + 0.5*dt*k1u, yp + 0.5*dt*k1v, tempo + dt*0.5);
      k2v = v(xp + 0.5*dt*k1u, yp + 0.5*dt*k1v, tempo + dt*0.5);
      k3u = u(xp + 0.5*dt*k2u, yp + 0.5*dt*k2v, tempo + dt*0.5);
      k3v = v(xp + 0.5*dt*k2u, yp + 0.5*dt*k2v, tempo + dt*0.5);
      k4u = u(xp + dt*k3u, yp + dt*k3v, tempo + dt);
      k4v = v(xp + dt*k3u, yp + dt*k3v, tempo + dt);
      // RK4
      (*itp)->set_particle_x(rk4(dt, xp, k1u, k2u, k3u, k4u));
      (*itp)->set_particle_y(rk4(dt, yp, k1v, k2v, k3v, k4v));	  

      
      start1 = clock();
      cparticle = M->search_particle_cell((*itp), number_of_levels - 1);
      if(cparticle == NULL){
	/*for(int k = number_of_levels - 1; k > 0; k--){
	  V = M->siblings((cparticle_bck));
	  if (V != NULL){
	    it = M->merge(V, cparticle_bck);
	    V->clear();
	  }
	  cparticle_bck = M->search_particle_cell(particle_bck, k-1);
	  }*/
	      
	for(int k = number_of_levels - 2; k >= 0; k--){
	  cparticle = M->search_particle_cell((*itp), k);
	  if(cparticle != NULL){
	     break;
	  }
	}
      }
      
      if(cparticle->get_cell_level() != number_of_levels - 1){
	//refina
	for(int k = cparticle->get_cell_level(); k < number_of_levels - 1; k++){
	  cparticle = M->search_particle_cell((*itp), k);
	  if(cparticle != NULL){
	    M->split(cparticle);
	  }
	}
      }
      //Neste ponto a particula esta no nivel mais fino
      (*itp)->set_particle_level(number_of_levels - 1);
      cparticle = M->search_particle_cell((*itp), (*itp)->get_particle_level());
      cparticle->set_cell_with_particle(1);
    }

    // Engrossamento so ocorre de number_of_levels != 1
    // percorre lista de celulas do nivel mais fino
    if(number_of_levels != 1){
      l = M->get_list_cell_by_level(number_of_levels - 1);
      it = l->begin();
      
      
      while (it != l->end()){
	V = M->siblings((*it));
	if (V != NULL && V->at(0)->get_cell_with_particle() == 0 && V->at(1)->get_cell_with_particle() == 0 && V->at(2)->get_cell_with_particle() == 0 && V->at(3)->get_cell_with_particle() == 0){
	  int xv = V->at(0)->get_cell_x()/2;
	  int yv = V->at(0)->get_cell_y()/2;
	  it = M->merge(V, (*it));
	  cparticle_bck = M->search(xv, yv, number_of_levels - 2);  
	  V->clear();
	  for(int k = number_of_levels - 2; k > 0; k--){
	    V = M->siblings((cparticle_bck));
	    if (V != NULL){
	      xv = V->at(0)->get_cell_x()/2;
	      yv = V->at(0)->get_cell_y()/2;
	      itv = M->merge(V, cparticle_bck);
	      cparticle_bck = M->search(xv, yv, k-1);  
	      V->clear();
	    }
	  }
	}
	else
	  it++;
      }
    }
    
    finish1 = clock();
    cputime += double(finish1 - start1)/CLOCKS_PER_SEC;
    //cout << ct << " CPU time = " << double(finish1 - start1)/CLOCKS_PER_SEC << endl;
    ct++;
    
    M->initialize_var(&u, &v, &phi, tempo, tempo);
    
    vmax = M->max_propriedades();
    dt = 0.5*h/max(vmax[0],vmax[1]);
    
    //if(ct%50 == 0)
    //M->print_silo(ct, &P);
    //M->get_hash_table()->print_information();
    //cout << dt << " " << tf << " " << tempo << " " << fabs(tf - tempo) << " " << tempo + dt << endl;
    dt = min(dt, fabs(tf - tempo));
    
    tempo += dt;    
    
  }
  finish = clock();
  //M->print_silo(ct, &P);
  printf("%E %E %E\n", double(finish0 - start0)/CLOCKS_PER_SEC, cputime, double(finish - start)/CLOCKS_PER_SEC);
  return 0;
}


