#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>

class PetscFachada {
public:
  PetscFacade();
  
  ~PetscFacade();
  
  Vec createVector(int size);
  
  Mat createMatrix(int rows, int cols); 
  
  void destroyVector(Vec& vec); 
  
  void destroyMatrix(Mat& mat);

  double norma_max(Vec& vec);
};



