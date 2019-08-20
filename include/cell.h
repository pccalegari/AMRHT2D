#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

class cell {
 private:
  int x, y;
  char kind [2]; // it could be bottom left (bl), bottom right (br), top left (tl), top right (tr)
  char last_kind [2]; // it could be bottom left (bl), bottom right (br), top left (tl), top right (tr) before last split
  int level;//0, 1, 2, 3, 4, ... from the courser to the finest  
  //double delta_x, delta_y;
  list <cell *>::iterator pointer_to_list;//ponteiro para a célula na lista de células da malha
  
 public:
  cell ();
  //initially this->kind and this->last_kind is equal to kind (= bl)
  cell (int x, int y, char kind [], int level);
  //constructor with parameter that specify the kind cell just before a split
  cell (int x, int y, char kind[], char last_kind[], int level);
  int get_cell_x();
  int get_cell_y();
  char * get_cell_kind();
  char * get_cell_last_kind();
  int get_cell_level();
  //double get_cell_delta_x();
  //double get_cell_delta_y();
  void set_cell_pointer_to_list(list<cell *>::iterator p);
  list<cell *>::iterator get_cell_pointer_to_list();
  cell ** split ();
  void print_cell ();
  void set_cell_kind(char kind[]);
  
};
