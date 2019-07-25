#include "hash_table.h"

hash_table::hash_table (int size){
  printf ("tamanho da tabela: %d\n\n", size);
  H = new vector<list <cell *> *>;
  list<cell *> * l;
  //(H->get_allocator()).allocate(size);
  for (int i = 0; i < size; i++){
    l = new list<cell *>;
    H->push_back(l);
  }
    
  //for (vector<list<cell *> *>::iterator it = H->begin(); it != H->end(); it++)
  //*it = new list<cell *>;
  number_cell = 0;
  load_factor = 0.;
}

void hash_table::insert(cell * c){
  int index = hash_function(c->get_cell_x(), c->get_cell_y(), c->get_cell_level());
  //printf ("index = %d\n", index);
  H->at(index)->push_back(c);
  number_cell++;
  load_factor = ((double) number_cell)/(H->size());
}

void hash_table::remove(cell * c){
  int index = hash_function(c->get_cell_x(), c->get_cell_y(), c->get_cell_level());
  for (list<cell *>::iterator it = H->at(index)->begin(); it != H->at(index)->end(); it++)
    if ((*it)->get_cell_x() == c->get_cell_x() && (*it)->get_cell_y() == c->get_cell_y() &&
	(*it)->get_cell_level() == c->get_cell_level()){
      H->at(index)->erase(it);
      break;
    }
  number_cell--;
  load_factor = ((double) number_cell)/(H->size());
}

cell * hash_table::search(int x, int y, int level){
  cell * founded_c = NULL;
  int index = hash_function(x, y, level);
  for (list<cell *>::iterator it = H->at(index)->begin(); it != H->at(index)->end(); it++)
    if ((*it)->get_cell_x() == x && (*it)->get_cell_y() == y && (*it)->get_cell_level() == level){
      founded_c = *it;
      break;
    }
  return founded_c;
}

int hash_table::hash_function(int x, int y, int level){
  //dimensão do domínio discretizado por nível h: n_h x n_h
  //n_0 x n_0 = 32 x 32
  //n_1 x n_1 = 64 x 64
  //n_2 x n_2 = 128 x 128
  //n_3 x n_3 = 256 x 256
  //n_4 x n_4 = 512 x 512
  int n[5] = {32, 64, 128, 256, 512};
  int key = x + y*n[level];
  for (int i = 0; i <= level - 1; i++)
    key += ((n[i] * n[i]) % H->size());//INI_SIZE_TABLE;
  int index = key % H->size();//INI_SIZE_TABLE;
  return index;
}

int hash_table::get_number_cell (){
  return number_cell;
}

double hash_table::get_load_factor(){
  return load_factor;
}

void hash_table::print_hash_table() {
  for (unsigned int i = 0; i < H->size();/*INI_SIZE_TABLE;*/ i++)
    if (!H->at(i)->empty()) {
      for (list<cell *>::iterator it = H->at(i)->begin(); it != H->at(i)->end(); it++) {
	/*assert((*(*it)->get_cell_pointer_to_list())->get_cell_x() == (*it)->get_cell_x() &&
	       (*(*it)->get_cell_pointer_to_list())->get_cell_y() == (*it)->get_cell_y() &&
	       (*(*it)->get_cell_pointer_to_list())->get_cell_level() == (*it)->get_cell_level());*/
	cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
      }
    }
}

unsigned int hash_table::get_hash_table_size() {
  return H->size();
}
