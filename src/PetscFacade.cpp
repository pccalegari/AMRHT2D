#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include "PetscFacade.h"

PetscFacade::PetscFacade() {
  PetscInitialize(nullptr, nullptr, nullptr, nullptr);
}

PetscFacade::~PetscFacade() {
  PetscFinalize();
}

Vec PetscFacade::createVector(int size) {
  Vec vec;
  VecCreate(PETSC_COMM_WORLD, &vec);
  VecSetSizes(vec, PETSC_DECIDE, size);
  VecSetFromOptions(vec);
  return vec;
}

Mat PetscFacade::createMatrix(int rows, int cols) {
  Mat mat;
  MatCreate(PETSC_COMM_WORLD, &mat);
  MatSetSizes(mat, PETSC_DECIDE, PETSC_DECIDE, rows, cols);
  MatSetFromOptions(mat);
  return mat;
}

void PetscFacade::destroyVector(Vec& vec) {
  VecDestroy(&vec);
}

void PetscFacade::destroyMatrix(Mat& mat) {
  MatDestroy(&mat);
}

double PetscFacade::norma_max(Vec& vec){
  double norm;
  VecNorm(vec, NORM_INFINITY, &norm);
  return norm;

}


