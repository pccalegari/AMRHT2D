#include "double_hash_table.h"

double_hash_table::double_hash_table (int size){
  //printf ("tamanho da tabela double HT: %d\n\n", size);
  H = new vector<list <double_cell *> *>;
  list<double_cell *> * l;
  //(H->get_allocator()).allocate(size);
  for (int i = 0; i < size; i++){
    l = new list<double_cell *>;
    H->push_back(l);
  }
    
  //for (vector<list<double_cell *> *>::iterator it = H->begin(); it != H->end(); it++)
  //*it = new list<double_cell *>;
  number_double_cell = 0;
  load_factor = 0.;
}

void double_hash_table::insert(double_cell * c){
  unsigned int index = double_hash_function(c->get_double_cell_x(), c->get_double_cell_y());
  //printf ("index = %d\n", index);
  H->at(index)->push_back(c);
  number_double_cell++;
  load_factor = ((double) number_double_cell)/(H->size());
}

void double_hash_table::remove(double_cell * c){
  unsigned int index = double_hash_function(c->get_double_cell_x(), c->get_double_cell_y());
  for (list<double_cell *>::iterator it = H->at(index)->begin(); it != H->at(index)->end(); it++)
    if (/*same key double cell*/(*it)->get_double_cell_x() == c->get_double_cell_x() && (*it)->get_double_cell_y() == c->get_double_cell_y()){
      H->at(index)->erase(it);
      break;
    }
  number_double_cell--;
  load_factor = ((double) number_double_cell)/(H->size());
}

double_cell * double_hash_table::search(double x, double y){
  double_cell * founded_c = NULL;
  unsigned int index = double_hash_function(x, y);
  for (list<double_cell *>::iterator it = H->at(index)->begin(); it != H->at(index)->end(); it++)
    if (/*same key double cell*/(*it)->get_double_cell_x() == x && (*it)->get_double_cell_y() == y){
      founded_c = *it;
      break;
    }
  return founded_c;
}

unsigned int double_hash_table::double_hash_function(double x, double y){
  double A = 0.6180339887;
  double B = 0.4237652378;
  double Ax, By;
  unsigned int m = H->size();
  Ax = A*x - floor(A*x);
  By = B*y - floor(B*y);
  unsigned int indexx, indexy;
  indexx = (unsigned int) (m * Ax);
  indexy = (unsigned int) (m * By);
  unsigned int index = (indexx + indexy) % m;
  return index;
}

int double_hash_table::get_number_double_cell (){
  return number_double_cell;
}

double double_hash_table::get_load_factor(){
  return load_factor;
}

void double_hash_table::print_double_hash_table() {
  for (unsigned int i = 0; i < H->size();/*INI_SIZE_TABLE;*/ i++)
    if (!H->at(i)->empty()) {
      for (list<double_cell *>::iterator it = H->at(i)->begin(); it != H->at(i)->end(); it++) {
	/*assert((*(*it)->get_double_cell_pointer_to_list())->get_double_cell_x() == (*it)->get_double_cell_x() &&
	       (*(*it)->get_double_cell_pointer_to_list())->get_double_cell_y() == (*it)->get_double_cell_y() &&
	       (*(*it)->get_double_cell_pointer_to_list())->get_double_cell_level() == (*it)->get_double_cell_level());*/
	cout << "(" << (*it)->get_double_cell_x() << ", " << (*it)->get_double_cell_y() << "): " << (*it)->get_double_cell_delta_x() << " " << (*it)->get_double_cell_delta_y() << endl;
      }
    }
}

unsigned int double_hash_table::get_double_hash_table_size() {
  return H->size();
}


void double_hash_table::clean(){
  vector <list<double_cell *> *>::iterator it = H->begin();
  while (it != H->end()) {
    if (!((*it)->empty())) 
      (*it)->clear();
    delete (*it);
    it++;
  }
  H->clear();
  delete H;
}


int double_hash_table::number_of_collision(){
  int count = 0;
  unsigned int max_size = 0;
  vector <list<double_cell *> *>::iterator it = H->begin();
  while (it != H->end()) {
    if (max_size < (*it)->size())
      max_size = (*it)->size();
    if ((*it)->size() > 1)
      count += ((*it)->size() - 1);
    it++;
  }
  cout << "Max size in one position: " << max_size << endl;
  return count;
}


void double_hash_table::print_information() {
  cout << "**********Secondary Hash Table (double keys)**********" << endl;
  cout << "Number of position on the hash table (m): "<< H->size() << endl;
  cout << "Number of elements stored (n): " << number_double_cell << endl;
  cout << "Load factor (n/m): " << load_factor << endl;
  cout << "Number of collisions: " << number_of_collision() << endl;
  cout << "**************************************************" << endl;
}
