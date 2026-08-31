#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h" 
#include <vector>

#define PI 3.1415926535897

using namespace std;

double f (double x, double y) { 
  return (cos(2.0 * PI * x) * sin(2.0 * PI * y));
}

double df (double x, double y) {
  return (-8.0 * PI * PI * f(x,y));
}

int main (){
  
  cell * c;
  //int level = 0;
  int number_of_levels = 1;
  int nxb = 4;
  int nyb = 4;
  vector<int> * max_dimension_by_level;
    
  //double delta_x, delta_y;
  dominio * D1, * D2;
  int nelem = 0;
  int nare = 0;
  int nvert = 0;
  int ncell = 0;
  double xbegin, ybegin, xend, yend, xbegin2, ybegin2, xend2, yend2;
  xbegin = ybegin = 0.;
  xend = 50.0;
  yend = 100.0;
  xbegin2 = 50.0;
  ybegin2 = 0.0;
  xend2 = yend2 = 100.0;
  D1 = new dominio (xbegin, ybegin, xend2, yend2);
  D2 = new dominio (xbegin2, ybegin2, xend2, yend2);
  mesh * M1, * M2;

  max_dimension_by_level = new vector<int>;
    
  //srand(time(NULL));
  srand(12345);
  
  /**********initiallize mesh*******************/
  for (int i = 0; i < number_of_levels; i++)
    max_dimension_by_level->push_back((nxb * pow(2,i)) * (nyb * pow(2,i)));
  M1 = new mesh(D1, number_of_levels, nxb ,nyb, max_dimension_by_level);
  //M2 = new mesh(D2, number_of_levels, nxb, nyb, max_dimension_by_level);
  
  /*********************************************/
  /* Nivel Base - Nivel 0 */
  for (int y = 0; y < nxb; y++)
    for (int x = 0; x < nyb; x++){
      c = new cell(x, y, 0,-1);
      M1->insert(c);
      //M2->insert(c);
    }
    
  /**************Refinement*************************/
  list <cell *> * l = M1->get_list_cell_by_level (0);
  list <cell *>::iterator it = l->begin();
  //hash_table * Hnew = M1->get_hash_table();

  /* Caso teste com refinamento estatico */
  /*while (it != l->end()) {
    if (((*it)->get_cell_x() >= 1 && (*it)->get_cell_x() < 2) && ((*it)->get_cell_y () >= 1 && (*it)->get_cell_y() <=2)) {
      c = M1->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
      assert (c != NULL);
      it = M1->split_and_insert(c);
    }
    else
      it++;
      }*/

  /******************************************/
  /* Apos a criacao da malha gerar os coeficientes da matriz uma lista de pesos para cada celula da malha */

  M1->get_hash_table()->print_information();
  //M2->get_hash_table()->print_information();
  
  //M->create_unstructured_mesh(&f, &df);

  M1->triangular_mesh_refined(&ncell, &nelem, &nare, &nvert);
  //cout << "P " << ncell << " " << nelem << " " << nare << " " << nvert << endl;
  //M2->triangular_mesh(&ncell, &nelem, &nare, &nvert);
  //cout << "Q " << ncell << " " << nelem << " " << nare << " " << nvert << endl;
  //M->mesh_adress(ncell, elem);
  	  
  return 0;
}
