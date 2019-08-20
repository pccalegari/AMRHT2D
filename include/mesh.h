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
  int nxb;
  int nyb;
  vector <int> * max_dimension_by_level;
 public:
  mesh();
  mesh(dominio * D, int n_levels, int nx, int ny, vector <int> * max_dimension_by_level);
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

  void left_neighbours_fd (cell * c, list <cell *> * nb);
  void right_neighbours_fd (cell * c, list <cell *> * nb);
  void up_neighbours_fd (cell * c, list <cell *> * nb);
  void down_neighbours_fd (cell * c, list <cell *> * nb);
  void insert_neighbour(cell * c, list <cell *> * nb);
  
  list <cell *> * neighbours (cell * c);
  list <cell *> * neighbours_fd (cell * c);
  
  void calculation_exact_function (double (* f) (double_cell * c));
  void calculation_function (double (* f) (double x, double y), vector <double> * fvalue);
  int counting_mesh_cell();
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
  
  void create_unstructured_mesh(double (* f) (double x, double y), double (* df)(double x, double y));
  dominio * get_dominio();
};
