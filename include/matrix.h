#include <iostream>
#include <list>
#include <vector>
#include "hash_table.h"
#include "double_hash_table.h"
#include "dominio.h"
#include "cell.h"
#include "mesh.h"
#include <list>

using namespace std;

//class matrix {
 //private:
 // vector <int> * ncell;
 //public:
 // matrix();
 // matrix(dominio * D, int n_levels, int nx, int ny, vector <int> * max_dimension_by_level);
 // void insert(cell *);
 // list <cell *>::iterator remove(cell *);
 // cell * search (int, int, int);
 // hash_table * get_hash_table();
 // list <cell *> * get_list_cell_by_level (int level);
 // list <cell *>::iterator split_and_insert(cell *);
  void left_neighbours_idx (cell * c, list <cell *> * nb);
  //void right_neighbours (cell * c, list <cell *> * nb);
  //void up_neighbours (cell * c, list <cell *> * nb);
  //void down_neighbours (cell * c, list <cell *> * nb);
  list <int *> * neighbours_idx (int * idx);
  //dominio * get_dominio();
//};
