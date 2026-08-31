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

double sol_u(double x, double y, double t) {
  //return (sin(2.0 * PI * x) * sin(2.0 * PI * y));
  return exp(0.5*x)*exp(0.5*y);
}

double f_rhs(double x, double y) {
  //return (8.0 * PI * PI * sin(2.0 * PI * x) * sin(2.0 * PI * y));
  return -0.5*exp(0.5*x)*exp(0.5*y);
}

//extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
//extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
//extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

static char help[] = "Solve a Poisson 2D in uniform mesh!";

int main (int argc, char **args){

  int number_of_levels = 2, nxb = 2, nyb;
  int ncell, idv, level;
  dominio *D;
  cell c;
  cell * ppc;

  double xbegin, ybegin, xend, yend;
  xbegin = ybegin = 0.0;
  xend = yend = 1.0;
  nyb = nxb;
  D = new dominio(xbegin, ybegin, xend, yend);

  mesh *M;
  M = new mesh(D, number_of_levels, nxb, nyb);
  double dx, dy, xd, yd;
  list <cell *> *l;
  list <cell *> *lv = nullptr;
  list <cell *>::iterator it, kt;
  
  KSP ksp;
  PetscReal norm, h, aux, h2;
  Vec uexa, rhs, r, u;
  Mat A;
  PetscInt i, j, iglobal, N;//, col0[3], col1[4], col2[5];
  PC pc;
  //PetscScalar w0[3], w1[4], w2[5];
  PetscMPIInt size;
  PetscLogDouble t1, t2;

  h = (xend - xbegin)/(double)nxb;
  h2 = h*h;

  /* Enumeração e lista de vizinhas */
  ncell = M->counting_mesh_cell();
  N = ncell;
  printf("%d\n", ncell);
  //M->print_mesh();
  //ncellv = M->neighbours_all_cell();

  M->create_unstructured_mesh(&sol_u, 0);
  //M->get_hash_table()->print_information();

  cout << "Teste com matriz - PETSc - Malha uniforme" << endl;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "Este eh um exemplo de um unico processador!");
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &N, NULL));
  
  /* Create right-size and solution vectors */
  PetscCall(VecCreate(PETSC_COMM_SELF, &u));
  PetscCall(VecSetSizes(u, PETSC_DECIDE, N));
  PetscCall(VecSetFromOptions(u));
  PetscCall(VecSet(u, 0.0));
  PetscCall(PetscObjectSetName((PetscObject)u, "Approx. Solution"));
  PetscCall(VecDuplicate(u, &rhs));
  PetscCall(PetscObjectSetName((PetscObject)rhs, "Right hand side"));
  PetscCall(VecDuplicate(rhs, &uexa));

  for (int kl = 0; kl < number_of_levels; kl++) {
    l = M->get_list_cell_by_level(kl);

    dx = fabs(xend - xbegin) / (nxb * pow(2, kl));
    dy = fabs(yend - ybegin) / (nyb * pow(2, kl));

    for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      level = (*it)->get_cell_level();
      xd = xbegin + ((i + 0.5) * dx);
      yd = ybegin + ((j + 0.5) * dy);

      iglobal = (*it)->get_cell_index();
      aux = (PetscScalar)(sol_u(xd, yd, 0));
      VecSetValues(uexa, 1, &iglobal, &aux, INSERT_VALUES);
      aux = h2 * (PetscScalar)(f_rhs(xd, yd));
      VecSetValues(rhs, 1, &iglobal, &aux, INSERT_VALUES);
      // Condição de Contorno Dirichlet - variável no centro
      if(i == 0){
        aux = 2.0 * (PetscScalar)(sol_u(xbegin, yd, 0));
        VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      }
      if(i == nxb-1){
        aux = 2.0 * (PetscScalar)(sol_u(xend, yd, 0));
        VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      }
      if(j == 0){
        aux = 2.0 * (PetscScalar)(sol_u(xd, ybegin, 0));
        VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      }
      if(j == nyb-1){
        aux = 2.0 * (PetscScalar)(sol_u(xd, yend, 0));
        VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
      }
    }
  }
   
  PetscCall(VecAssemblyBegin(uexa));
  PetscCall(VecAssemblyEnd(uexa));
  PetscCall(VecAssemblyBegin(rhs));
  PetscCall(VecAssemblyEnd(rhs));

  /* Create Matrix */
  PetscCall(MatCreate(PETSC_COMM_SELF, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  vector <PetscInt> idx;
  vector <PetscScalar> wp;

  for (int kl = 0; kl < number_of_levels; kl++) {
    l = M->get_list_cell_by_level(kl);

    dx = fabs(xend - xbegin) / (nxb * pow(2, kl));
    dy = fabs(yend - ybegin) / (nyb * pow(2, kl));

    for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      level = (*it)->get_cell_level();
      xd = xbegin + ((i + 0.5) * dx);
      yd = ybegin + ((j + 0.5) * dy);

      iglobal = (*it)->get_cell_index();
      ppc = M->search(i, j, level);
      //list <cell *> * lnew;
      if (ppc != NULL) {
        lv = M->neighbours(ppc);
        //printf("%d ",iglobal);
        if(lv->size() == 2){
          wp.push_back(6.0);
        }
        else if(lv->size() == 3){
          wp.push_back(5.0);
        }
        else{
          wp.push_back(4.0);
        }
        idx.push_back(iglobal);

        for (list <cell *>::iterator ite = lv->begin(); ite != lv->end(); ite++) {
          idv = (*ite)->get_cell_index();
          idx.push_back(idv);
          wp.push_back(-1.0);
          //printf("%d ", idv);
        }
        //std::cout << "Tamanho: " << idx.size() << " " << wp.size() << "\n";
        //delete lnew;
        //std::cout << "Vazio?: " << lv->size() << "\n";
        //printf("\n");
      }
      else{
        cout << "There is no such cell!" << endl;
      }
      lv->clear();
      PetscCall(MatSetValues(A, 1, &iglobal, idx.size(), idx.data(), wp.data(), INSERT_VALUES));
      idx.clear();
      idx.shrink_to_fit();
      wp.clear();
      wp.shrink_to_fit();
      //cout << "Tamanho: " << idx.size() << " " << wp.size() << "\n";
    }
  }


   PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
   PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
   PetscCall(MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));

   //PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
   //PetscCall(VecView(rhs, PETSC_VIEWER_STDOUT_WORLD));

   /*Set Solver linear*/

   PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
   PetscCall(KSPSetOperators(ksp, A, A));
   PetscCall(KSPGetPC(ksp, &pc));
   PetscCall(KSPSetType(ksp, KSPBCGS));

   PetscCall(PCSetType(pc, PCILU));
   //PetscCall(KSPSetType(ksp, KSPPREONLY));
   PetscCall(KSPSetTolerances(ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
   PetscCall(KSPSetFromOptions(ksp));

   /* Solver */
   PetscCall(PetscTime(&t1));
   PetscCall(KSPSolve(ksp, rhs, u));
   PetscCall(PetscTime(&t2));
   PetscPrintf(PETSC_COMM_WORLD, "Tempo de execução: %f segundos\n", t2 - t1);

   PetscCall(KSPGetSolution(ksp, &u));
   PetscCall(KSPGetRhs(ksp, &rhs));
   PetscCall(VecDuplicate(rhs, &r));
   PetscCall(KSPGetOperators(ksp, &A, NULL));
   
   PetscCall(MatMult(A, u, r));
   PetscCall(VecAXPY(r, -1.0, rhs));
   PetscCall(VecNorm(r, NORM_2, &norm));
   PetscCall(PetscPrintf(PETSC_COMM_SELF, "Residual norm %g\n", (double)norm));
   //for(i = 0; i < N*N; i++){
   //printf("%f %f\n", u[i], uexa[i]);
   //  }

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


