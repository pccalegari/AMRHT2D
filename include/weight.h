#include <iostream>
#include <list>
#include <vector>
//#include "cell.h"

using namespace std;

/*
class weight {
 private:
  int cell_index;
  double w;
      
 public:
  weight();
  weight(int index, double w);
  int get_weight_index();
  double get_weight_w();
  void set_weight_index(int index);
  void set_weight_w(double w);
  void add_weight_w(double w);
  void left_neighbours_idx (cell *, vector <int > *);
  //list <int >  neighbours_idx (int * idx);
};
*/

//  int get_weight_index();
//  double get_weight_w();
//  void set_weight_index(int index);
//  void set_weight_w(double w);
//  void add_weight_w(double w);
  void left_neighbours_idx (cell * c, vector <int > * idx);
  vector <int> neighbours_idx (cell * c, vector <int> * idx);
