#include<iostream>
#include "hash_table.h"
#include "math.h"

hash_table::hash_table (int size, int max_level, int n_colunas_max_level){
  //printf ("tamanho da tabela: %d\n", size);
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

  this->n_colunas_max_level = n_colunas_max_level;
  this->max_level = max_level;
  
}

void hash_table::insert(cell * c){
  unsigned int index = hash_function(c->get_cell_x(), c->get_cell_y(), c->get_cell_level());
  assert (index < H->size());
  //printf ("index = %d\n", index);
  H->at(index)->push_back(c);
  number_cell++;
  load_factor = ((double) number_cell)/(H->size());
  //cout << "enumera = " << number_cell << " " << c->get_cell_x() << " " << c->get_cell_y() << endl;
}

void hash_table::remove(cell * c){
  unsigned int index = hash_function(c->get_cell_x(), c->get_cell_y(), c->get_cell_level());
  assert (index < H->size());
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
  unsigned int index = hash_function(x, y, level);
  assert(index < H->size());
  for (list<cell *>::iterator it = H->at(index)->begin(); it != H->at(index)->end(); it++)
    if ((*it)->get_cell_x() == x && (*it)->get_cell_y() == y && (*it)->get_cell_level() == level){
      founded_c = *it;
      break;
    }
  return founded_c;
}

/*unsigned int hash_table::hash_function_old(int x, int y, int level){
  int n[6] = {32, 64, 128, 256, 512, 1024};
  int key = x + y*n[level];
  for (int i = 0; i <= level - 1; i++)
    key += ((n[i] * n[i]) % H->size());
  int index = key % H->size();
  return index;
  }*/

unsigned int hash_table::hash_function(int x, int y, int level){

  int index, key, i, j;
  int diff_level = this->max_level - level;
  int pow_2_diff_level = pow(2, diff_level);
  double ra, kr;
  
  ra = (sqrt(5) - 1.0)/2.0;
  i = x * pow_2_diff_level; //encontra a abscissa correspondente a x no nível mais fino
  j = y * pow_2_diff_level; //encontra a ordenada correspondente a y no nível mais fino
  key = (this->n_colunas_max_level * j + i); //encontra a numeração da célula (i, j, l_max)
                                             //considerando numeração da esquerda para a direira e
                                             //de cima para baixo (e começando do 0)
  /* Metodo da Divisao */
  index = key % H->size();
  
  /* Metodo de Fibonacci */
  kr = key*ra;
  index = (int)(H->size()*(kr - (int)(kr)));
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
	cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << ", " << (*it)->get_cell_index() << endl;
      }
    }
  cout << "Número maximo de níveis da malha: " << max_level << endl;
  cout << "Número de colunas do nível mais fino: " << n_colunas_max_level << endl;
							       
}

unsigned int hash_table::get_hash_table_size() {
  return H->size();
}

int hash_table::number_of_collision(){
  int count = 0, ct;
  unsigned int max_size = 0;
  
  vector <list<cell *> *>::iterator it = H->begin();
  while (it != H->end()) {
    //cout << "# cell position " << (*it)->size() << endl;
    if (max_size < (*it)->size())
      max_size = (*it)->size();
    if ((*it)->size() > 1) {
      count += ((*it)->size() - 1);
    }
    it++;
  }
  for(unsigned int k = 0; k <= max_size; k++){
    it = H->begin();
    ct = 0;
    while (it != H->end()) {
      if((*it)->size() == k)
	ct++;
      it++;
    }
    cout << k << " " << ct << endl;
  }
  cout << "Max size in one position: " << max_size << endl;
  
  return count;
}

void hash_table::print_information() {
  cout << "**********Main Hash Table (integer keys)**********" << endl;
  cout << "Number of position on the hash table (m): "<< H->size() << endl;
  cout << "Number of elements stored (n): " << number_cell << endl;
  cout << "Load factor (n/m): " << load_factor << endl;
  cout << "Number of collisions: " << number_of_collision() << endl;
  cout << "**************************************************" << endl;
}
