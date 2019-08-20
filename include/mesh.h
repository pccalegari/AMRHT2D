#include <iostream>
#include <list>
#include <vector>
#include "silo.h"
#include "hash_table.h"
#include "double_hash_table.h"
#include "dominio.h"

using namespace std;

class mesh {
 private:
  dominio * D;
  hash_table * H;
  vector<list <cell *> *> * l;
  int number_of_levels;
  vector <int> * max_dimension_by_level;
 public:
  mesh();
  mesh(dominio * D, int n_levels, vector <int> * max_dimension_by_level);
  void insert(cell *);
  list <cell *>::iterator remove(cell *);
  cell * search (int, int, int);
  void print_mesh();
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
  
  void calculation_exact_function (double (* f) (double_cell * c));
  void create_unstructured_mesh(double (* f) (double x, double y, double t), double tempo);
  dominio * get_dominio();

  //Recebe uma célula c e devolve uma lista de células irmãs de c somente se existirem. Caso contrário, devolve NULL. As células irmãs de c são c mais outras três células produzidas com o refinamento de uma célula.
  vector <cell *> * siblings (cell * c);

  //Recebe uma lista de células irmãs e agrupa todas elas. Recebe também o ponteiro para a célula na lista de células que está sendo percorrida (importante para atualizar iterador para a função que chama merge). Devolve um iterador para a próxima célula que não é nenhuma daquelas que foram agrupadas.
  list <cell *>::iterator merge (vector <cell *> * V, cell * c);
 
  
};
