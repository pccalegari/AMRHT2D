#include "petsc.h"

//#undef __funct__
//#define __funct__ "main"

int main(int argc, char **args){
  PetscErrorCode ierr;
  PetscMPIInt rank;

  PetscInitialize(&argc, &args, (char *)0, PETSC_NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  ierr = PetscPrintf(PETSC_COMM_SELF, "Hello by procs %d!\n", rank);CHKERRQ(ierr);
  ierr=PetscFinalize();
  return 0;
}
