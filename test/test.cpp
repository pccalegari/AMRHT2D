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
  int number_of_levels = 2;
  int nxb = 4;
  int nyb = 4;
  vector<int> * max_dimension_by_level;
  vector<double> *rhs, *u, *ff, *A;
  vector<int> *JA, *IA;
  
  //double delta_x, delta_y;
  dominio * D;
  int ncell = 0;
  int ncellv = 0;
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.;
  xend = yend = 1.0;
  D = new dominio (xbegin, ybegin, xend, yend);
  mesh * M;

  //delta_x = delta_y = 0.05;
  max_dimension_by_level = new vector<int>;
  rhs = new vector<double>;
  u = new vector<double>;
  ff = new vector<double>;
    
  //srand(time(NULL));
  srand(12345);
  
  /**********initiallize mesh*******************/
  for (int i = 0; i < number_of_levels; i++)
    max_dimension_by_level->push_back((nxb * pow(2,i)) * (nyb * pow(2,i)));
  M = new mesh(D, number_of_levels, nxb ,nyb, max_dimension_by_level);
  /*********************************************/
  /* Nivel Base - Nivel 0 */
  for (int y = 0; y < nxb; y++)
    for (int x = 0; x < nyb; x++){
      c = new cell(x, y, 0,-1);
      M->insert(c);
    }
    
  /**************Refinement*************************/
  list <cell *> * l = M->get_list_cell_by_level (0);
  list <cell *>::iterator it = l->begin();
  //hash_table * Hnew = M->get_hash_table();

  /* Caso teste com refinamento estatico */
  while (it != l->end()) {
    if (((*it)->get_cell_x() >= 1 && (*it)->get_cell_x() < 2) && ((*it)->get_cell_y () >= 1 && (*it)->get_cell_y() <=2)) {
      c = M->search((*it)->get_cell_x(), (*it)->get_cell_y(), (*it)->get_cell_level());
      assert (c != NULL);
      it = M->split_and_insert(c);
    }
    else
      it++;
  }

  /******************************************/
  /* Apos a criacao da malha gerar os coeficientes da matriz uma lista de pesos para cada celula da malha */

  M->get_hash_table()->print_information();
  
  M->create_unstructured_mesh(&f, &df);

  ncell = M->counting_mesh_cell();
  ncellv = M->neighbours_all_cell();
  cout << ncell << " " << ncellv << endl;
  vector <cell *> * elem;
  elem = new vector <cell *> (ncell);
  IA = new vector <int> (ncell+1,-1);
  JA = new vector <int> (ncellv,-1);
  A = new vector <double> (ncellv,0.0);
  
  M->mesh_adress(ncell, elem);
  M->create_matrix_df (A, JA, IA, ncell, ncellv, elem);
  M->print_matrix(IA, JA, A, ncell, ncellv);
  
  M->calculation_function (&df, ff);
  M->calculation_function (&f, u);
  M->rhs_dirichlet_boundary_conditions (&df, rhs);
  
  delete rhs;
  delete u;
  delete ff;
  delete elem;
  delete IA;
  delete A;
  delete JA;
	  
  return 0;
}
