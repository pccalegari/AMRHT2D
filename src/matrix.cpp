#include "mesh.h"
#include <vector>
#include <list>

void matrix::left_neighbours_idx (cell * c, list <int *> * idx) {
    int x = c->get_cell_x();
    int y = c->get_cell_y();
    int level = c->get_cell_level();
    int xl, yl, k;

    if (x != 0) {
        int l;
    if (x % 2 == 0){
        for (l = level - 1; l >= 0; l--){
            xl = (int) (x / pow (2, level - l));
            yl = (int) (y / pow (2, level - l));
            xl = xl - 1;
            //printf ("(%d, %d): %d)\n", xl, yl, i);
            cell * found = H->search (xl, yl, l);
            if (found != NULL) {
                k = (*found)->get_cell_index();
                idx->insert (idx->begin(), k);
                break;
            }
        }
    }
    if (x % 2 == 1 || (x % 2 == 0 && l == -1)) {//não encontrou uma célula à esquerda mais grossa
        int xx;
        for (l = level; l < number_of_levels; l++) {
            xx = (int) (pow (2, l - level) * x);
            for (xl = xx - 1; xl < xx; xl++){
                for (int j = 0; j <= (int) pow (2, l - level) - 1; j++) {
                    yl = (int) (pow (2, l - level) * y) + j;
                    //printf ("(%d, %d): %d)\n", xl, yl, i);
                    cell * found = H->search (xl, yl, l);
                    if (found != NULL) {
                        k = (*found)->get_cell_index();
                        idx->insert (idx->begin(), k);
                    }
                }
            }
        }
    }
}
}
/*
void mesh::right_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;

  if (x != nxb*pow(2,level)) {
    int i;
    if (x % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	xl = xl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);
	  //cout << "FOUND in coarser" << endl;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (x % 2 == 0 || (x % 2 == 1 && i == -1)) {//não encontrou uma célula à esquerda mais grossa
      for (i = level; i < number_of_levels; i++) {
	xl = (int) (pow (2, i - level) * x);
	xl = xl + pow(2, i - level); //procuramos pelos vizinhos à direita de c e a referência de c é o canto inferior esquerdo. Temos que "pular" as células nível atual (+ fino) dentro do próprio c
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  yl = (int) (pow (2, i - level) * y) + j;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->begin(), found);
	    //cout << "FOUND in finer" << endl;
	    //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  }
	}
      }
    }
  }
}

void mesh::down_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;

  if (y != 0) {
    int i;
    if (y % 2 == 0){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl - 1;
	//printf ("(%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL){
	  nb->insert (nb->begin(), found);
	  //found = H->search(xl-1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("DL.: (%d, %d): %d)\n", xl-1, yl, i);
	  //}
	  //found = H->search(xl+1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("DR.: (%d, %d): %d)\n", xl+1, yl, i);
	  //}
	  //cout << "FOUND in coarser" << endl;
	  //printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (y % 2 == 1 || (y % 2 == 0 && i == -1)) {//não encontrou uma célula abaixo mais grossa
      for (i = level; i < number_of_levels; i++) {
	yl = (int) (pow (2, i - level) * y);
	yl = yl - 1;
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  xl = (int) (pow (2, i - level) * x) + j;
	  //printf ("(%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->begin(), found);
	    //cout << "FOUND in finer" << endl;
	    //printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	    //found = H->search (xl-1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf ("DL.: (%d, %d): %d)\n", xl-1, yl, i);
	    //}
	    //found = H->search (xl+1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf ("DR.: (%d, %d): %d)\n", xl+1, yl, i);
	  }

	}
      }
    }
  }

}

void mesh::up_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;

  if (y != nyb*pow(2,level)) {
    int i;
    if (y % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);

	  //found = H->search(xl-1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("UL.: (%d, %d): %d)\n", xl-1, yl, i);
	  //}
	  //found = H->search(xl+1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("UR.: (%d, %d): %d)\n", xl+1, yl, i);
	  //}
	  //cout << "FOUND in coarser" << endl;
	  //printf ("U.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (y % 2 == 0 || (y % 2 == 1 && i == -1)) {//não encontrou uma célula à esquerda mais grossa
      for (i = level; i < number_of_levels; i++) {
	yl = (int) (pow (2, i - level) * y);
	yl = yl + pow(2, i - level); //procuramos pelos vizinhos à direita de c e a referência de c é o canto inferior esquerdo. Temos que "pular" as células nível atual (+ fino) dentro do próprio c
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  xl = (int) (pow (2, i - level) * x) + j;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->begin(), found);
	    //found = H->search (xl-1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf("UL.: (%d, %d): %d)\n", xl-1, yl, i);
	    //}
	    //found = H->search (xl+1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf("UR.: (%d, %d): %d)\n", xl+1, yl, i);
	    //}
	    //cout << "FOUND in finer" << endl;
	    //printf ("U.: (%d, %d): %d)\n", xl, yl, i);
	  }
	}
      }
    }
  }

}
*/

list <cell *> * matrix::neighbours_idx (cell * c) {
  list <int *> * idx = new list <int*>;
  //nb->insert(nb->begin(), c);
  this->left_neighbours_idx(c, idx);
  //this->right_neighbours(c, nb);
  //this->down_neighbours(c, nb);
  //this->up_neighbours(c, nb);

  return idx;
}
