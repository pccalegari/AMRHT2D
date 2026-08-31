#include <iostream>
#include <list>
#include <vector>
#include "hash_table.h"
#include "double_hash_table.h"
#include "dominio.h"
#include "cell.h"

using namespace std;

class matrix {
 private:
  vector <int> * ncell + 1;
 public:
  matrix();
  matrix(dominio * D, int n_levels, int nx, int ny, vector <int> * max_dimension_by_level);
  void insert(cell *);
  list <cell *>::iterator remove(cell *);
  cell * search (int, int, int);
  hash_table * get_hash_table();
  list <cell *> * get_list_cell_by_level (int level);
  list <cell *>::iterator split_and_insert(cell *);

  /*Recebe uma célula c e o nível mais fino que c (ou mesmo nível de c) de seus vizinho(s) à esquerda.
    Devolve uma lista com seu(s) vizinho(s). Essa lista NÃO PODE SER NULA.*/
  list <cell *> * left_finest_or_same_level_neighbours (cell * c, int level_k);

  void left_neighbours (cell * c, list <cell *> * nb);
  void right_neighbours (cell * c, list <cell *> * nb);
  void up_neighbours (cell * c, list <cell *> * nb);
  void down_neighbours (cell * c, list <cell *> * nb);
  list <cell *> * neighbours (cell * c);
  dominio * get_dominio();
};
