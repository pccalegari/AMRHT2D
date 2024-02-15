#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "mesh.h" 
#include <vector>

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

#define PI 3.1415926535897

using namespace std;

double sol_u(double x, double y) { 
  return (sin(2.0 * PI * x) * sin(2.0 * PI * y));
}

double f_rhs(double x, double y) {
  return (8.0 * PI * PI * sin(2.0 * PI * x) * sin(2.0 * PI * y));
}

//extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
//extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
//extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

static char help[] = "Solve a Poisson 2D!";

int main (int argc, char **args){

  int number_of_levels = 2, nxb = 4, nyb = 4, N = nxb, index;

  dominio *D;
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.0;
  xend = yend = 1.0;
  D = new dominio(xbegin, ybegin, xend, yend);

  mesh *M;
  M = new mesh(D, number_of_levels, nxb, nyb);
  double dx, dy, xd, yd;
  list <cell *> *l;
  list <cell *>::iterator it;
  
  KSP ksp;
  PetscReal norm, x, y, aux;
  Vec uexa, u, rhs, r;
  Mat A;
  PetscInt i, j, iglobal, col0[3], col1[4], col2[5];
  PC pc;
  PetscScalar w0[3], w1[4], w2[5];
  PetscMPIInt size = 1;

  /* Geração da malha */
  for(int i = 0; i < number_of_levels - 1; i++){
    l = M->get_list_cell_by_level(i);
    dx = fabs(xend - xbegin)/(nxb * pow(2, i));
    dy = fabs(yend - ybegin)/(nyb * pow(2, i));

    it = l->begin();

    while(it != l->end()){
      xd = xbegin + (((*it)->get_cell_x()) * dx);
      yd = ybegin + (((*it)->get_cell_y()) * dy);
      if (f_rhs(xd+0.5*dx, yd+0.5*dy) < 0.999 - i*0.05 && f_rhs(xd+0.5*dx, yd+0.5*dy) > - 0.999 + i*0.05){
      //if (df(xd, yd, 0.0) < 0.9 && df(xd, yd, 0.0) > - 0.9){
      	it = M->split(*it);
      }
      else
	it++;
      
    }
   
  }
    
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "Este eh um exemplo de um unico processador!");
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &nxb, NULL));
  
  
  /* Create right-size and solution vectors */
  //u = Petsc.createVector(N*N);
  PetscCall(VecCreate(PETSC_COMM_SELF, &u));
  PetscCall(VecSetSizes(u, PETSC_DECIDE, nxb*nxb));
  PetscCall(VecSetFromOptions(u));
  PetscCall(PetscObjectSetName((PetscObject)u, "Approx. Solution"));
  //rhs = Petsc.createVector(N*N);
  PetscCall(VecDuplicate(u, &rhs));
  PetscCall(PetscObjectSetName((PetscObject)rhs, "Right hand side"));
  PetscCall(VecDuplicate(rhs, &uexa));
  PetscCall(VecSet(u, 0.0));

  dx = fabs(xend - xbegin) / nxb;
  dy = fabs(yend - ybegin) / nyb;
  index = 0;
  for (int i = 0; i < number_of_levels; i++) {
    
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    
    //it = l->begin();
    for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
      xd = xbegin + (((*it)->get_cell_x() + 0.5) * dx);
      yd = ybegin + (((*it)->get_cell_y() + 0.5) * dy);
      (*it)->set_cell_index(index);
      printf("%d\n", (*it)->get_cell_index()); 

      index++;
      
    }
  }
  
  for(j = 0; j < nyb; j++){
    y = ybegin + (j + 0.5)*dx;
    for(i = 0; i < nxb; i++){
      x = xbegin + (i + 0.5)*dy;
      /*printf("%f %f", x, y);*/
      iglobal = i + j*nxb;
      aux = (PetscScalar)(sol_u(x, y));
      VecSetValues(uexa, 1, &iglobal, &aux, INSERT_VALUES);
      // Alterar para incluir contorno nao nulo 
      aux = dx*dy * (PetscScalar)(f_rhs(x, y));
      VecSetValues(rhs, 1, &iglobal, &aux, INSERT_VALUES);
    }
  }
  
  /* Alteracao do lado direito por causa da variavel no centro da celula */
     
   for(j = 0; j < nyb; j++){
      y = ybegin + (j + 0.5)*dx;
      // i = 0
      aux = 2.0 * (PetscScalar)(sol_u(xbegin, y));
      iglobal = j*nxb;
      VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      // i = N-1
      aux = 2.0 * (PetscScalar)(sol_u(xend, y));
      iglobal = nxb - 1 + j*nxb;
      VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      // j = 0 
      aux = 2.0 * (PetscScalar)(sol_u(y, ybegin));
      iglobal = j;
      VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      // j = N-1
      aux = 2.0 * (PetscScalar)(sol_u(y, yend));
      iglobal = j + (nyb-1)*nxb;
      VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      
   }
 
   PetscCall(VecAssemblyBegin(uexa));
   PetscCall(VecAssemblyEnd(uexa));
   PetscCall(VecAssemblyBegin(rhs));
   PetscCall(VecAssemblyEnd(rhs));

   /* Create Matrix */
  
   w0[0] = 6.0;
   w0[1] = -1.0;
   w0[2] = -1.0;
    
   w1[0] = 5.0;
   w1[1] = -1.0;
   w1[2] = -1.0;
   w1[3] = -1.0;

   w2[0] = 4.0;
   w2[1] = -1.0;
   w2[2] = -1.0;
   w2[3] = -1.0;
   w2[4] = -1.0;

   PetscCall(MatCreate(PETSC_COMM_SELF, &A));
   PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N*N, N*N));
   PetscCall(MatSetFromOptions(A));
   PetscCall(MatSetUp(A));
   
   for (j = 1; j < N - 1; j++) {
      for (i = 1; i < N - 1; i++){
         iglobal = i + j*N; 
         col2[0] = i + j*N;
         col2[1] = i + (j+1)*N;
         col2[2] = i + (j-1)*N;
         col2[3] = i - 1 + j*N;
         col2[4] = i + 1 + j*N;
         PetscCall(MatSetValues(A, 1, &iglobal, 5, col2, w2, INSERT_VALUES));
      }
   }
   //linha 0 canto do domínio
   iglobal = 0;
   col0[0] = 0;
   col0[1] = N;
   col0[2] = 1;
   PetscCall(MatSetValues(A, 1, &iglobal, 3, col0, w0, INSERT_VALUES));
   // linha N-1 canto do domínio
   iglobal = N-1;
   col0[0] = N-1;
   col0[1] = N-2;
   col0[2] = 2*N-1;
   PetscCall(MatSetValues(A, 1, &iglobal, 3, col0, w0, INSERT_VALUES));
   //linha N*N-N canto do domínio
   iglobal = N*(N-1);
   col0[0] = N*(N-1);
   col0[1] = 1 + N*(N-1);
   col0[2] = N*(N-2);
   PetscCall(MatSetValues(A, 1, &iglobal, 3, col0, w0, INSERT_VALUES));
   // linha N-1 canto do domínio
   iglobal = N*N-1;
   col0[0] = N*N-1;
   col0[1] = N*N-2;
   col0[2] = N-1+(N-2)*N;
   PetscCall(MatSetValues(A, 1, &iglobal, 3, col0, w0, INSERT_VALUES));
  
   for(j = 1; j < N -1; j++){
      iglobal = j*N; 
      col1[0] = j*N;
      col1[1] = (j+1)*N;
      col1[2] = (j-1)*N;
      col1[3] = 1 + j*N;
      PetscCall(MatSetValues(A, 1, &iglobal, 4, col1, w1, INSERT_VALUES));
      iglobal = N-1 + j*N; 
      col1[0] = N-1 + j*N;
      col1[1] = N-1 + (j+1)*N;
      col1[2] = N-1 + (j-1)*N;
      col1[3] = N-2 + j*N;
      PetscCall(MatSetValues(A, 1, &iglobal, 4, col1, w1, INSERT_VALUES));
      iglobal = j; 
      col1[0] = j;
      col1[1] = j + N;
      col1[2] = j - 1;
      col1[3] = j + 1;
      PetscCall(MatSetValues(A, 1, &iglobal, 4, col1, w1, INSERT_VALUES));
      iglobal = j+(N-1)*N; 
      col1[0] = j+(N-1)*N;
      col1[1] = j-1+N*(N-1);
      col1[2] = j+1+N*(N-1);
      col1[3] = j+(N-2)*N;
      PetscCall(MatSetValues(A, 1, &iglobal, 4, col1, w1, INSERT_VALUES));
   }
      
   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
   PetscCall(MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));

   /*PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
   PetscCall(VecView(rhs, PETSC_VIEWER_STDOUT_WORLD));*/


   /*Set Solver linear*/
   
   PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
   PetscCall(KSPSetOperators(ksp, A, A));
   PetscCall(KSPGetPC(ksp, &pc));
   PetscCall(PCSetType(pc, PCJACOBI));
   PetscCall(KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
   PetscCall(KSPSetFromOptions(ksp));

   /* Solver */
   
   PetscCall(KSPSolve(ksp, rhs, u));
   PetscCall(KSPGetSolution(ksp, &u));
   PetscCall(KSPGetRhs(ksp, &rhs));
   PetscCall(VecDuplicate(rhs, &r));
   PetscCall(KSPGetOperators(ksp, &A, NULL));
   
   PetscCall(MatMult(A, u, r));
   PetscCall(VecAXPY(r, -1.0, rhs));
   PetscCall(VecNorm(r, NORM_2, &norm));
   PetscCall(PetscPrintf(PETSC_COMM_SELF, "Residual norm %g\n", (double)norm));
   /*for(i = 0; i < N*N; i++){
     printf("%f %f\n", u[i], uexa[i]);
     }*/
   
   PetscCall(VecAXPY(u, -1.0, uexa));
   PetscCall(VecNorm(u, NORM_INFINITY, &norm));
   PetscCall(PetscPrintf(PETSC_COMM_SELF, "Erro norm %g\n", (double)norm));

   PetscCall(KSPDestroy(&ksp));
  
   PetscCall(VecDestroy(&r));
   
   PetscCall(VecDestroy(&rhs));
   PetscCall(MatDestroy(&A));
   PetscCall(VecDestroy(&u));
   PetscCall(VecDestroy(&uexa));
   PetscFunctionReturn(PETSC_SUCCESS);
   PetscCall(PetscFinalize());
   return 0;
 }


