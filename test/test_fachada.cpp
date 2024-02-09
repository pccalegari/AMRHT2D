#include "PetscFacade.h"

int main(int argc, char **argv) {
    PetscFacade petsc;

    // Uso básico
    Vec vector = petsc.createVector(10);
    Mat matrix = petsc.createMatrix(10, 10);

    double norma = petsc.norma_max(vector);
    printf("%f", norma);
    petsc.destroyVector(vector);
    petsc.destroyMatrix(matrix);

    return 0;
}
