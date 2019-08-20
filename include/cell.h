#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

class cell {
 private:
  int x, y;
  int level;//0, 1, 2, 3, 4, ... from the courser to the finest
  int index;
  //double delta_x, delta_y;
  list <cell *>::iterator pointer_to_list;//ponteiro para a célula na lista de células da malha
  
 public:
  cell ();
  //cell (int x, int y, int level, double delta_x, double delta_y);
  cell (int x, int y, int level, int index);
  int get_cell_x();
  int get_cell_y();
  int get_cell_level();
  int get_cell_index();
  void put_cell_index(int ivalue);
  //double get_cell_delta_x();
  //double get_cell_delta_y();
  void set_cell_pointer_to_list(list<cell *>::iterator p);
  list<cell *>::iterator get_cell_pointer_to_list();
  cell ** split ();
  void print_cell ();
  
};
