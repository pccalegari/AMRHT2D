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

  int number_of_levels = 1, nxb = 64, nyb = 64;
  int ncell;
  dominio *D;
  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.0;
  xend = yend = 1.0;
  D = new dominio(xbegin, ybegin, xend, yend);
  D->set_tbc_left(1);
  D->set_tbc_right(1);
  D->set_tbc_up(1);
  D->set_tbc_down(1);
  
  mesh *M;
  M = new mesh(D, number_of_levels, nxb, nyb);
  double dx, dy, xd, yd;
  list <cell *> *l, *lnb;
  list <cell *>::iterator it, kt;
  
  KSP ksp;
  PetscReal norm, aux;
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

  /* Enumeração e lista de vizinhas */
  ncell = M->counting_mesh_cell();
  printf("%d\n", ncell);

  for (i = 0; i < number_of_levels; i++) {
    l = M->get_list_cell_by_level(i);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    
    for (list <cell *>::reverse_iterator it = l->rbegin(); it != l->rend(); it++) {
      xd = xbegin + (((*it)->get_cell_x() + 0.5) * dx);
      yd = ybegin + (((*it)->get_cell_y() + 0.5) * dy);
      iglobal = (*it)->get_cell_index();
      aux = (PetscScalar)(sol_u(xd, yd));
      VecSetValues(uexa, 1, &iglobal, &aux, INSERT_VALUES);
      aux = dx*dy * (PetscScalar)(f_rhs(xd, yd));
      VecSetValues(rhs, 1, &iglobal, &aux, INSERT_VALUES);
      (*it)->set_cell_nb(M->neighbours(*it));
      //printf("Celula:\n");
      //(*it)->print_cell();
      //printf("Vizinhas:\n");
      //M->print_list((*it)->get_cell_nb());
    }
  }

  /* Alteracao do lado direito por causa da variavel no centro da celula */
  list <cell *> * bcs = M->boundary_cells_south(D->get_tbc_down());
  for(list <cell *>::iterator it = bcs->begin(); it != bcs->end(); it++){
    int l = (*it)->get_cell_level();
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    xd = xbegin + ((*it)->get_cell_x() + 0.5)*dx;
    aux = 2.0 * (PetscScalar)(sol_u(xd, ybegin));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bcn = M->boundary_cells_north(D->get_tbc_up());
  for(list <cell *>::iterator it = bcn->begin(); it != bcn->end(); it++){
    int l = (*it)->get_cell_level();
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    xd = xbegin + ((*it)->get_cell_x() + 0.5)*dx;
    aux = 2.0 * (PetscScalar)(sol_u(xd, yend));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bce = M->boundary_cells_east(D->get_tbc_left());
  for(list <cell *>::iterator it = bce->begin(); it != bce->end(); it++){
    int l = (*it)->get_cell_level();
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    yd = ybegin + ((*it)->get_cell_y() + 0.5)*dy;
    aux = 2.0 * (PetscScalar)(sol_u(xbegin, yd));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bcw = M->boundary_cells_west(D->get_tbc_right());
  for(list <cell *>::iterator it = bcw->begin(); it != bcw->end(); it++){
    int l = (*it)->get_cell_level();
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    yd = ybegin + ((*it)->get_cell_y() + 0.5)*dy;
    aux = 2.0 * (PetscScalar)(sol_u(xend, yd));
    iglobal = (*it)->get_cell_index();
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
   PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, ncell, ncell));
   PetscCall(MatSetFromOptions(A));
   PetscCall(MatSetUp(A));
   
   for (i = 0; i < number_of_levels; i++) {
     l = M->get_list_cell_by_level(i);
     
     for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
       iglobal = (*it)->get_cell_index();
       lnb = (*it)->get_cell_nb();
       if((*it)->get_cell_bc() == 0){
	 j = 0;
	 for(kt = lnb->begin(); kt != lnb->end(); kt++){
	   col2[j] = (*kt)->get_cell_index();
	   //printf("%d %d %d\n", iglobal, j, col2[j]);
	   j++;
	 }
	 PetscCall(MatSetValues(A, 1, &iglobal, j, col2, w2, INSERT_VALUES));
       }
       else{
	 j = 0;
	 int tam = lnb->size();
	 if(tam == 3){
	   for(kt = lnb->begin(); kt != lnb->end(); kt++){
	     col0[j] = (*kt)->get_cell_index();
	     j++;
	   }
	   PetscCall(MatSetValues(A, 1, &iglobal, j, col0, w0, INSERT_VALUES));
	 }
	 else{
	   for(kt = lnb->begin(); kt != lnb->end(); kt++){
	     col1[j] = (*kt)->get_cell_index();
	     j++;
	   }
	   PetscCall(MatSetValues(A, 1, &iglobal, j, col1, w1, INSERT_VALUES));
	 }
       }
     }
   }
   
   
   //linha 0 canto do domínio
   /*
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
   */
      
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


