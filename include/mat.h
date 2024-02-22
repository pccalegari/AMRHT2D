#include <iostream>
#include <list>
#include <vector>
#include "mesh.h"

using namespace std;

class mat {
   private:
      
 public:
  mat();
  mat(mesh * M);
  mesh * get_mesh();
};
