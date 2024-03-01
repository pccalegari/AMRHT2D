#include <iostream>
#include <list>
#include <vector>

using namespace std;

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
};
