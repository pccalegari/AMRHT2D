#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h"
#include "particle.h"

#define PI 3.1415926535897

using namespace std;

double f (double x, double y) { 
  return cos(2 * PI * x) * sin(2 * PI * y);
}

double df (double x, double y, double tempo) {
  return -4 * PI * PI * f(tempo * x, tempo * y);
}

int main (){
  
  cell * c;

  list <particle *> P;

  particle * p = new particle();

  P.push_back(p);

  particle * q = P.back();

  q->print_particle();
  
  //int level = 0;
  int number_of_levels = 4;
  vector<int> * max_dimension_by_level;
  //double delta_x, delta_y;
  dominio * D;
  
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = -1.;
  xend = yend = 1.;
  D = new dominio (xbegin, ybegin, xend, yend);
  
  mesh * M;
  
  //delta_x = delta_y = 0.05;
  max_dimension_by_level = new vector<int>;
  
  //srand(time(NULL));
  srand(12345);
  
  /**********initiallize mesh*******************/
  for (int i = 0; i < number_of_levels; i++)
    max_dimension_by_level->push_back((32 * pow(2,i)) * (32 * pow(2,i)));
  M = new mesh(D, number_of_levels, max_dimension_by_level);
  /*********************************************/

  char kind [2];
  kind[0] = 'b'; kind[1] = 'l';
  for (int y = 0; y < 32; y++)
    for (int x = 0; x < 32; x++){
      c = new cell(x, y, kind, 0);
      M->insert(c);
    }
    
  /**************Refinement*************************/
  //list <cell *> * l = M->get_list_cell_by_level (0);
  //list <cell *>::iterator it = l->begin();
  //int r;
  //while (it != l->end()) {
  //  if (/*1*/((*it)->get_cell_x() >= 10 && (*it)->get_cell_x () <= 21)/* r = rand() % 3 == 0 || r % 2 == 1*/) {
  //   c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
  //    assert (c != NULL);
  //    it = M->split_and_insert(c);
  //  }
  //  else
  //    it++;
  //}

  //l = M->get_list_cell_by_level (1);
  //it = l->begin();
  //while (it != l->end()) {
  //  if (/*1*//*r = rand() % 3 == 0 ||*/ rand() % 3 == 1) {
  //    c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
  //    assert (c != NULL);
  //    it = M->split_and_insert(c);
  //  }
  //  else
  //    it++;
  //}

  //l = M->get_list_cell_by_level (2);
  //it = l->begin();
  //while (it != l->end()) {
  //  if (/*1*//*r = rand() % 3 == 0 ||*/ rand() % 3 == 1) {
  //    c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
  //    assert (c != NULL);
  //    it = M->split_and_insert(c);
  //  }
  //  else
  //    it++;
  //}

  //l = M->get_list_cell_by_level (3);
  //it = l->begin();
  //while (it != l->end()) {
  //  if (/*1*//*rand() % 3 == 0 ||*/ rand() % 3 == 1) {
  //    c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
  //    assert (c != NULL);
  //    it = M->split_and_insert(c);
  //  }
  //  else
  //    it++;
  //}
  
  //l = M->get_list_cell_by_level (4);
  //it = l->begin();
  //while (it != l->end()) {
  //  if (/*1*//*rand() % 3 == 0 ||*/ rand() % 3 == 1) {
  //    c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
  //    assert (c != NULL);
  //    it = M->split_and_insert(c);
  //  }
  //  else
  //    it++;
  //}
  /******************************************/
  
  M->get_hash_table()->print_information();

  
  double dx, dy, xd, yd;
  list <cell *> * l;
  list <cell *>::iterator it;
  double tempo;

  xbegin = M->get_dominio()->get_xbegin();
  ybegin = M->get_dominio()->get_ybegin();
  xend = M->get_dominio()->get_xend();
  yend = M->get_dominio()->get_yend();
  
  for (tempo = -1; tempo <= 0; tempo += 0.005){
    //printf ("%lf\n", tempo);
    //refinement & merge
    for (int i = number_of_levels - 1; i >= 0 ; i--) {
      l = M->get_list_cell_by_level(i);
      
      if (l->size() > 0){
	it = l->begin();
	//(*it)->print_cell();
	//cout << endl;
	dx = fabs(xend - xbegin) / (32 * pow(2, i));
	dy = fabs(yend - ybegin) / (32 * pow(2, i));
	//int contador = 0;
	while (it != l->end()){
	  
	  xd = xbegin + ((*it)->get_cell_x() * dx);//(delta_0 / pow(2, i)));
	  yd = ybegin + ((*it)->get_cell_y() * dy);//(delta_0 / pow(2, i)));

	  assert(((*it)->get_cell_kind()[0] == 'b' || (*it)->get_cell_kind()[0] == 't') && ((*it)->get_cell_kind()[1] == 'l' || (*it)->get_cell_kind()[1] == 'r'));

	  //refinement
	  if ((df(xd + (dx / 2.), yd + (dy / 2.), tempo) <= 25 && df(xd + (dx / 2.), yd + (dy / 2.), tempo) >= 18) ||
	      (df(xd + (dx / 2.), yd + (dy / 2.), tempo) >= -25 && df(xd + (dx / 2.), yd + (dy / 2.), tempo) <= -18)){
	    if (i < number_of_levels - 1 /*se i é diferente do último nível, então pode refinar*/){
	      //c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
	      //assert(c != NULL);
	      //assert(c == *it);
	      it = M->split_and_insert(*it);
	    }
	    else {
	      it++;
	    }
	    //***
	  }
	  else{
	    //it++; //código sem merge
	    assert(i == (*it)->get_cell_level());
	    //Mais uma implementação do merge ***************************
	    if (i > 0) {//se i é diferente do primeiro nível, então pode fazer merge
	      vector <cell *> * V;
	      V = M->siblings((*it));
	      if (V != NULL){
		//cout << V->size() << endl;
		it = M->merge(V, *it);
		V->clear();
	      }
	      else it++;
	    }
	    else it++;
	  }
	}
      }
    }
    M->create_unstructured_mesh(&df, tempo);
    M->get_hash_table()->print_information();
  }
  
  M->create_unstructured_mesh(&df, tempo);
  M->get_hash_table()->print_information();
  
  
  cell * cnew = new cell (47, 16, kind, 2);
  
  hash_table * Hnew = M->get_hash_table();

  //Hnew->print_hash_table();

  if (Hnew->search(47, 16, 2) != NULL) {

    list <cell *> * lnew = M->neighbours (cnew);
    
    for (list <cell *>::iterator it = lnew->begin(); it != lnew->end(); it++) {
      (*it)->print_cell();
    }
  }
  else
    {
      cout << "There is no such cell!" << endl;
    }
  return 0;
}
