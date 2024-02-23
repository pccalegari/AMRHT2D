#include "weight.h"

weight::weight(){
  cell_index = 0;
  w = 0.0;
}

weight::weight(int index, double w){
  this->cell_index = index;
  this->w = w;
}

int weight::get_weight_index(){
  return cell_index;
}

double weight::get_weight_w(){
  return w;
}

void weight::set_weight_index(int index){
  cell_index = index;
}

void weight::set_weight_w(double wp){
  w += wp;
}


  
