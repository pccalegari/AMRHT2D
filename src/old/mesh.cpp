#include <cmath>
#include <cassert>
#include "mesh.h"

mesh::mesh(){
  l = NULL;
  number_of_levels = -1;
  max_dimension_by_level = NULL;
}

mesh::mesh(int n_levels, vector <int> * max_dimension_by_level){
  //o tamanho da tabela hash vai depender do máximo número de células do nível mais fino, isto é:
  //(1/100)*max_dim_by_level[n_level - 1] => 10 % de max_dim_by_level[n_level - 1]
  H = new hash_table((int)(0.1 * max_dimension_by_level->at(n_levels - 1)));
  l = new vector<list<cell *> * >;
  for (int i = 0; i < n_levels; i++){
    list<cell *> * l_by_level= new list <cell *>;
    l->push_back(l_by_level); 
  }
  number_of_levels = n_levels;
  this->max_dimension_by_level = max_dimension_by_level;
}

void mesh::insert(cell * c){
  c->set_cell_pointer_to_list((l->at(c->get_cell_level())->insert(l->at(c->get_cell_level())->begin(), c)));
  H->insert(c); 
}

void mesh::remove(cell * c) {
  l->at(c->get_cell_level())->erase(c->get_cell_pointer_to_list());
  H->remove(c);
  delete c;
}

cell * mesh::search(int x, int y, int level) {
  return H->search(x, y, level);
}

void mesh::print_mesh() {
  for (int i = 0; i < number_of_levels; i++)
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++)
      cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
  
}

hash_table * mesh::get_hash_table() {
  return H;
}


void mesh::split_and_insert(cell * c) {
  cell ** V;
  V = c->split();
  remove(c);
  insert(V[0]);
  insert(V[1]);
  insert(V[2]);
  insert(V[3]);
  free (V);
}

list <cell *> * mesh::left_finest_or_same_level_neighbours(cell * c, int level_k) {
  list <cell*> * l = new list <cell *>;
  //O valor em level_k DEVE SER maior ou igual que c->get_cell_level();
  int k = level_k - c->get_cell_level();
  assert (k >= 0);
  int number_neighbours_by_x = (int) pow(2, k);
  int number_neighbours_by_y = (int) pow(2, k);
  int _2_to_k = number_neighbours_by_x;
  int i, j, x, y;
  x = c->get_cell_x();
  y = c->get_cell_y();
  cell * neighbour_c;
  for (j = 0; j < number_neighbours_by_y; j++) {
    for (i = 0; i < number_neighbours_by_x; i++) {
      //(x - 1) determina os vizinhos à esquerda de c
      neighbour_c = search (_2_to_k * (x - 1) + i, _2_to_k * y + j, level_k);
      assert(neighbour_c != NULL);
      l->push_back(neighbour_c);
    }
  }
  return l;
}

void mesh::create_unstructured_mesh() {
  //4 pontos para cada celula (haverá pontos para duas ou mais células, com isso, o fator de carga de internal_HT será menor que 1
  //int x, y;
  double xd, yd;
  //por enquanto as variáveis abaixo (a, b, delta) descrevem o domínio. (a,b) é o ponto onde começa o domínio e delta é o fator de discretização no nível mais grosso (nível 0)
  double a, b, delta;
  a = b = 0;
  delta = 1.;
  double_hash_table * internal_HT = new double_hash_table(4 * H->get_number_cell());
  for (int i = 0; i < number_of_levels; i++) {
    printf ("Nível: %d\n", i);
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      //cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
      xd = a + ((*it)->get_cell_x() * (delta / pow(2, i)));
      yd = b + ((*it)->get_cell_y() * (delta / pow(2, i)));
      
    }
  
  } 
}
