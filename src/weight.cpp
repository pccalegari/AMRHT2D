#include "cell.h"
#include "weight.h"

/*
weight::weight(){
  cell_index = 0;
  w = 0.0;
}

weight::weight(int index, double w){
  this->cell_index = index;
  this->w = w;
}

*/
/*
int get_weight_index(){
  return cell_index;
}

double get_weight_w(){
  return w;
}
*/

//void set_weight_index(int index){
//  cell_index = index;
//}

//void set_weight_w(double wp){
//  w = wp;
//}

//void add_weight_w(double wp){
//  w += wp;
//}

void left_neighbours_idx (cell * c, vector <int> * idx) {
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
                    idx->push_back(k);//->insert (idx->begin(), k);
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
                            idx->push_back(k);//->insert (idx->begin(), k);
                        }
                    }
                }
            }
        }
    }

}

vector <int> neighbours_idx (cell * c, vector <int> * idx) {
  idx = new vector <int>;
  //nb->insert(nb->begin(), c);
  this->left_neighbours_idx(c, idx);
  //this->right_neighbours(c, nb);
  //this->down_neighbours(c, nb);
  //this->up_neighbours(c, nb);

  return idx;
}
*/
  
