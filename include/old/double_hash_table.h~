#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include "cell.h"

//#define INI_SIZE_TABLE 1000000

using namespace std;

class hash_table {
 private:
  vector <list<cell *> *> * H;
  //<list<cell *> *> * H;
  int number_cell;
  double load_factor;
  
 public:
  hash_table (int size);
  void insert(cell *);
  void remove(cell *);
  //Recebe: a chave de uma célula, ou seja, a sua posição (x, y) e o seu nível h
  //Devolve: ou NULL (se a célula não existe); ou o endereço da célula (se a célula existe)
  cell * search(int, int, int);
  //Recebe: a chave de uma célula, ou seja, a sua posição (x, y) e o seu nível h
  //Devolve: um índice na tabela hash H
  int hash_function(int, int, int);
  int get_number_cell ();
  double get_load_factor();
  unsigned int get_hash_table_size();
  void print_hash_table();
};
