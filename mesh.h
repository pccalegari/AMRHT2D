#include <iostream>
#include <list>
#include <vector>
#include "silo.h"
#include "hash_table.h"
#include "double_hash_table.h"
#include "dominio.h"
#include "particle.h"
//#include "weight.h"

using namespace std;

class mesh {
 private:
  dominio * D;
  hash_table * H;
  vector<list <cell *> *> * l;
  int number_of_levels;
  int nxb;
  int nyb;
  vector <int> * max_dimension_by_level;
 public:
  mesh();
  //mesh(dominio * D, int n_levels, int nx, int ny, vector <int> * max_dimension_by_level);
  mesh(dominio * D, int n_levels, int nx, int ny);
  void insert(cell *);
  list <cell *>::iterator remove(cell *);
  cell * search (int, int, int);
  cell * search_particle_cell (particle *, int);
  void print_mesh();
  hash_table * get_hash_table();
  int get_nxb();
  int get_nyb();
  list <cell *> * get_list_cell_by_level (int level);
  list <cell *>::iterator split(cell *);
  cell ** split_return_new_cells(cell * c);

  /*Recebe uma célula c e o nível mais fino que c (ou mesmo nível de c) de seus vizinho(s) à esquerda.
    Devolve uma lista com seu(s) vizinho(s). Essa lista NÃO PODE SER NULA.*/
  list <cell *> * left_finest_or_same_level_neighbours (cell * c, int level_k);

  int left_neighbours (cell * c, list <cell *> * nb);
  int right_neighbours (cell * c, list <cell *> * nb);
  int up_neighbours (cell * c, list <cell *> * nb);
  int down_neighbours (cell * c, list <cell *> * nb);

  void left_neighbours_fd (cell * c, list <cell *> * nb);
  void right_neighbours_fd (cell * c, list <cell *> * nb);
  void up_neighbours_fd (cell * c, list <cell *> * nb);
  void down_neighbours_fd (cell * c, list <cell *> * nb);
  void insert_neighbour(cell * c, list <cell *> * nb);
  
  list <cell *> * neighbours (cell * c);
  list <cell *> * neighbours_fd (cell * c);

  list <cell *> * neighbours_e(cell * c);
  list <cell *> * neighbours_w(cell * c);
  list <cell *> * neighbours_s(cell * c);
  list <cell *> * neighbours_n(cell * c);
    
  void calculation_exact_function (double (* f) (double_cell * c));
  void calculation_function (double (* f) (double x, double y), vector <double> * fvalue);
  int counting_mesh_cell();
  
  list <cell *> * boundary_cells_south(int tbc);
  list <cell *> * boundary_cells_north(int tbc);
  list <cell *> * boundary_cells_east(int tbc);
  list <cell *> * boundary_cells_west(int tbc);
  void print_cell_list(list <cell *> * bc);
  void print_weight_list(list <weight> bc);
  
  void triangular_mesh(int * ncellp, int * nelemp, int * narep, int * nvertp);
  void triangular_mesh_refined(int * ncellp, int * nelemp, int * narep, int * nvertp);
  int search_left_neighbours(cell * c);
  int search_down_neighbours(cell * c);
  int search_right_neighbours(cell * c);
  int search_up_neighbours(cell * c);
  int search_vertice(int i, int j, int n, int id, int l, vector <cell *> * vert);
  void ordena_vertice(int nvert, vector <int> * vind, vector <float> * vx, vector <float> * vy);
  int search_vertice_index(float x, float y, int nvert, vector <float> * vx, vector <float> * vy, vector <int> * vind);
  int neighbours_all_cell();
  void mesh_adress(int ncell, vector<cell *> * elem);
  void create_matrix_df (vector<double> * IA, vector<int> * JA, vector<int> * II, int ncell, int ncellv, vector <cell*> * elem);
  double normmax(vector <double> * x, int ncell);
  double erron2(vector <double> * x, vector <double> * xa, int ncell);
  void gs_csr(vector <double> * x, vector <double> * A, vector <int> * IA, vector <int> * JA, vector <double> * b, int ncell);
  void print_matrix(vector <int> * IA, vector <int> * JA, vector <double> * A, int ncell, int ncellv);
  void rhs_dirichlet_boundary_conditions (double (*f) (double x, double y), vector<double> * fvalue);
  
  void create_unstructured_mesh(double (* f) (double x, double y, double t), double tempo);
  void print_silo(int ct, list <particle *> *P);
  void initialize_var(double (* u) (double x, double y, double t), double (* v) (double x, double y, double t), double (*phi) (double x, double y, double t), double tempo, double tempo0);
  vector <double> max_propriedades();
  //Recebe uma célula c e devolve uma lista de células irmãs de c somente se existirem. Caso contrário, devolve NULL. As células irmãs de c são c mais outras três células produzidas com o refinamento de uma célula.
  vector <cell *> * siblings (cell * c);

  //Recebe uma lista de células irmãs e agrupa todas elas. Recebe também o ponteiro para a célula na lista de células que está sendo percorrida (importante para atualizar iterador para a função que chama merge). Devolve um iterador para a próxima célula que não é nenhuma daquelas que foram agrupadas.
  list <cell *>::iterator merge (vector <cell *> * V, cell * c);

  dominio * get_dominio();
};
