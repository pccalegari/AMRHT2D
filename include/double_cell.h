#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

class double_cell {
 private:
  double x, y;
  double delta_x, delta_y;
  double variable;
  list <double_cell *>::iterator pointer_to_list;//ponteiro para a célula na lista de células da malha
  int pointer_to_vector;
  
 public:
  double_cell ();
  //double_cell (int x, int y, int level, double delta_x, double delta_y);
  double_cell (double x, double y, double dx, double dy);
  double get_double_cell_x();
  double get_double_cell_y();
  double get_double_cell_delta_x();
  double get_double_cell_delta_y();
  //double get_double_cell_delta_x();
  //double get_double_cell_delta_y();
  void set_double_cell_pointer_to_list(list<double_cell *>::iterator p);
  list<double_cell *>::iterator get_double_cell_pointer_to_list();
  double_cell ** split ();
  void print_double_cell ();
  void set_double_cell_pointer_to_vector(int index);
  int get_double_cell_pointer_to_vector();
  void set_center_variable (double (* f) (double x, double y));
};
