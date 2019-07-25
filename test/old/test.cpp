#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include "mesh.h" 

using namespace std;

int main (){

  cell * c;
  int level = 0;
  int number_of_levels = 4;
  vector<int> * max_dimension_by_level;
  //double delta_x, delta_y;
  mesh * M;

  //delta_x = delta_y = 0.05;
  max_dimension_by_level = new vector<int>;
  
  
  /**********initiallize mesh*******************/
  for (int i = 0; i < number_of_levels; i++)
    max_dimension_by_level->push_back((32 * pow(2,i)) * (32 * pow(2,i)));
  M = new mesh(number_of_levels, max_dimension_by_level);
  /*********************************************/

  for (int y = 0; y < 32; y++)
    for (int x = 0; x < 32; x++){
      c = new cell(x, y, level);
      M->insert(c);
    }

  //level++;
  //delta_x = delta_y = delta_x/2.;
  
  for (int y = 0; y < 32; y++)
    for (int x = 14; x <= 17; x++){
      c = M->search(x, y, level);
      assert (c != NULL);
      M->split_and_insert(c);
    }

  list <cell *> * l_neighbours = M->left_finest_or_same_level_neighbours(M->search(18,14,0), 1);
  assert (l_neighbours->size() == 4);
  for (list <cell *>::iterator it = l_neighbours->begin(); it != l_neighbours->end(); it++) {
    (*it)->print_cell();
  }
  delete l_neighbours;
    

  /*assert(M->search(0,2,0) != NULL);
  printf ("(%d, %d):%d\n", M->search(31,31,0)->get_cell_x(), M->search(31,31,0)->get_cell_y(), M->search(31,31,0)->get_cell_level());

  printf ("(%d, %d)\n\n", (*(M->search(31,31,0)->get_cell_pointer_to_list()))->get_cell_x(), (*(M->search(31,31,0)->get_cell_pointer_to_list()))->get_cell_y());

  assert(M->search(0,1024, 5) == NULL);

  M->get_hash_table()->print_hash_table();
  
  M->print_mesh();*/

  printf ("número de elementos = %d\nfator de carga = %f\n",M->get_hash_table()->get_number_cell(), M->get_hash_table()->get_load_factor());

  DBfile *dbfile = NULL;
  /* Open the Silo file */
  dbfile = DBCreate("data/basic.silo", DB_CLOBBER, DB_LOCAL,"Comment about the data", DB_HDF5);
  if(dbfile == NULL)
    {
      fprintf(stderr, "Could not create Silo file!\n");
      return -1;
    }

  

  /* Add other Silo calls here. */
  /* Close the Silo file. */
  DBClose(dbfile);  

  M->create_unstructured_mesh();
    
  return 0;
}
