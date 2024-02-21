#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include "weight.h"

using namespace std;

class cell {
 private:
  int x, y;
  int level;//0, 1, 2, 3, 4, ... from the courser to the finest
  int index;
  int bc;
  double phi, phi0, velu, velv;
  list <weight> * wpoisson;
  list <cell *> * nb;
  int cp;
  list <cell *>::iterator pointer_to_list;//ponteiro para a célula na lista de células da malha
  
 public:
  cell ();
  cell (int x, int y, int level);
  int get_cell_x();
  int get_cell_y();
  int get_cell_level();
  int get_cell_index();
  int get_cell_bc();
  int get_cell_with_particle();
  void set_cell_index(int ivalue);
  void set_cell_bc(int tbc);
  double get_cell_velu();
  double get_cell_velv();
  double get_cell_phi();
  double get_cell_phi0();
  list <cell *> * get_cell_nb();
  list <weight> * get_cell_wpoisson();
  void set_cell_nb(list <cell *> * vnb);
  void set_cell_wpoisson(list <weight> * lw);
  void set_cell_velu(double uvalue);
  void set_cell_velv(double vvalue);
  void set_cell_phi(double phivalue);
  void set_cell_phi0(double phi0value);
  void set_cell_with_particle(int wparticle);
  void set_cell_pointer_to_list(list<cell *>::iterator p);
  list<cell *>::iterator get_cell_pointer_to_list();
  cell ** split ();
  void print_cell ();
};
