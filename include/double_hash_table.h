#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <cmath>
#include "double_cell.h"

//#define INI_SIZE_TABLE 1000000

using namespace std;

class double_hash_table {
 private:
  vector <list<double_cell *> *> * H;
  //<list<double_cell *> *> * H;
  int number_double_cell;
  double load_factor;
  
 public:
  double_hash_table (int size);
  void insert(double_cell *);
  void remove(double_cell *);
  //Recebe: a chave de uma célula, ou seja, a sua posição (x, y) e o seu nível h
  //Devolve: ou NULL (se a célula não existe); ou o endereço da célula (se a célula existe)
  double_cell * search(double, double);
  //Recebe: a chave de uma célula, ou seja, a sua posição (x, y) e o seu nível h
  //Devolve: um índice na tabela hash H
  unsigned int double_hash_function(double, double);
  int get_number_double_cell ();
  double get_load_factor();
  unsigned int get_double_hash_table_size();
  void print_double_hash_table();
  void clean();
  int number_of_collision();
  void print_information();
};
