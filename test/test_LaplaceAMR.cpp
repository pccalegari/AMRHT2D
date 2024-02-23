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

double Dif(double x, double y){
  return 1.0;
}

double f_rhs(double x, double y) {
  return (8.0 * PI * PI * sin(2.0 * PI * x) * sin(2.0 * PI * y));
}

//extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
//extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
//extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

static char help[] = "Solve a Poisson 2D!";

int main (int argc, char **args){

  int number_of_levels = 1, nxb = 32, nyb = 64;
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
  double dx, dy, xd, yd, dx2, dy2, De, Dw, Ds, Dn, Dp;
  list <cell *> *l;
  list <weight> MP;
  list <weight>::iterator wt;
  list <cell *> * lvw, * lve, * lvn, * lvs;
  list <cell *>::iterator it, kt;
  
  KSP ksp;
  PetscReal norm, aux;
  Vec uexa, u, rhs, r;
  Mat A;
  PetscInt i, j, k, iglobal, col2[5];
  PC pc;
  PetscScalar w2[5];
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
  
  
  

  /* Enumeração e lista de vizinhas */
  ncell = M->counting_mesh_cell();
  printf("%d\n", ncell);

  /* Create right-size and solution vectors */
  //u = Petsc.createVector(N*N);
  PetscCall(VecCreate(PETSC_COMM_SELF, &u));
  PetscCall(VecSetSizes(u, PETSC_DECIDE, ncell));
  PetscCall(VecSetFromOptions(u));
  PetscCall(PetscObjectSetName((PetscObject)u, "Approx. Solution"));
  //rhs = Petsc.createVector(N*N);
  PetscCall(VecDuplicate(u, &rhs));
  PetscCall(PetscObjectSetName((PetscObject)rhs, "Right hand side"));
  PetscCall(VecDuplicate(rhs, &uexa));
  PetscCall(VecSet(u, 0.0));
  
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
      aux = (PetscScalar)(f_rhs(xd, yd));
      VecSetValues(rhs, 1, &iglobal, &aux, INSERT_VALUES);
      (*it)->set_cell_lbw(M->neighbours_w(*it));
      (*it)->set_cell_lbe(M->neighbours_e(*it));
      (*it)->set_cell_lbn(M->neighbours_n(*it));
      (*it)->set_cell_lbs(M->neighbours_s(*it));
      //printf("Celula:\n");
      //(*it)->print_cell();
      //printf("Vizinhas:\n");
      //M->print_cell_list((*it)->get_cell_lbw());
      //M->print_cell_list((*it)->get_cell_lbs());
      //M->print_cell_list((*it)->get_cell_lbe());
      //M->print_cell_list((*it)->get_cell_lbn());
    }
  }


  /* Alteracao do lado direito por causa da variavel no centro da celula */
  /* Coeficientes da matriz associados as fronteiras */
  
  list <cell *> * bcs = M->boundary_cells_south(D->get_tbc_down());
  for(list <cell *>::iterator it = bcs->begin(); it != bcs->end(); it++){
    int l = (*it)->get_cell_level();
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    dy2 = dy*dy;
    xd = xbegin + ((*it)->get_cell_x() + 0.5)*dx;
    Ds = -Dif(xd, ybegin)/dy2;
    aux = Ds * 2.0 * (PetscScalar)(sol_u(xd, ybegin));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bcn = M->boundary_cells_north(D->get_tbc_up());
  for(list <cell *>::iterator it = bcn->begin(); it != bcn->end(); it++){
    int l = (*it)->get_cell_level();
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    dy2 = dy*dy;
    xd = xbegin + ((*it)->get_cell_x() + 0.5)*dx;
    Dn = -Dif(xd, yend)/dy2;
    aux = Dn * 2.0 * (PetscScalar)(sol_u(xd, yend));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bce = M->boundary_cells_east(D->get_tbc_right());
  for(list <cell *>::iterator it = bce->begin(); it != bce->end(); it++){
    int l = (*it)->get_cell_level();
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    dx2 = dx*dx;
    yd = ybegin + ((*it)->get_cell_y() + 0.5)*dy;
    De = -Dif(xend, yd)/dx2;
    aux = De * 2.0 * (PetscScalar)(sol_u(xend, yd));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }
  
  list <cell *> * bcw = M->boundary_cells_west(D->get_tbc_left());
  for(list <cell *>::iterator it = bcw->begin(); it != bcw->end(); it++){
    int l = (*it)->get_cell_level();
    dy = fabs(yend - ybegin) / (nyb * pow(2, l));
    dx = fabs(xend - xbegin) / (nxb * pow(2, l));
    dx2 = dx*dx;
    yd = ybegin + ((*it)->get_cell_y() + 0.5)*dy;
    Dw = -Dif(xbegin, yd)/dx2;
    aux = Dw * 2.0 * (PetscScalar)(sol_u(xbegin, yd));
    iglobal = (*it)->get_cell_index();
    VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
  }

  /* Cálculo dos pesos Poisson */
  
  for (k = 0; k < number_of_levels; k++) {
    l = M->get_list_cell_by_level(k);
    
    dx = fabs(xend - xbegin) / (nxb * pow(2, k));
    dy = fabs(yend - ybegin) / (nyb * pow(2, k));
    dx2 = dx*dx;
    dy2 = dy*dy;
    
    for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
      weight wp, we, ww, ws, wn;
      list <weight> lw;
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      xd = xbegin + ((i + 0.5) * dx);
      yd = ybegin + ((j + 0.5) * dy);
      De = -Dif(xd + 0.5*dx, yd)/dx2;
      Dw = -Dif(xd - 0.5*dx, yd)/dx2;
      Ds = -Dif(xd, yd - 0.5*dy)/dy2;
      Dn = -Dif(xd, yd + 0.5*dy)/dy2;
      lvw = (*it)->get_cell_lbw();
      lvs = (*it)->get_cell_lbs();
      lve = (*it)->get_cell_lbe();
      lvn = (*it)->get_cell_lbn();
            
      Dp = -(De + Dw + Ds + Dn);
           
      if(lvw->size() != 0){
	for(list <cell *>::iterator wt = lvw->begin(); wt != lvw->end(); wt++){
	  ww.set_weight_index((*wt)->get_cell_index());
	  ww.set_weight_w(Dw); //Precisa corrigir Dw varia no AMR
	  lw.insert(lw.end(), ww);
	}
      }
      else{
	Dp -= Dw;
      }
      if(lvs->size() != 0){
	for(list <cell *>::iterator wt = lvs->begin(); wt != lvs->end(); wt++){
	  ws.set_weight_index((*wt)->get_cell_index());
	  ws.set_weight_w(Ds); //Precisa corrigir Dw varia no AMR
	  lw.insert(lw.end(), ws);
	}
      }
      else{
	Dp -= Ds;
      }
      if(lve->size() != 0){
	for(list <cell *>::iterator wt = lve->begin(); wt != lve->end(); wt++){
	  we.set_weight_index((*wt)->get_cell_index());
	  we.set_weight_w(De); //Precisa corrigir Dw varia no AMR
	  lw.insert(lw.end(), we);
	  
	}
      }
      else{
	Dp -= De;
      }
      if(lvn->size() != 0){
	for(list <cell *>::iterator wt = lvn->begin(); wt != lvn->end(); wt++){
	  wn.set_weight_index((*wt)->get_cell_index());
	  wn.set_weight_w(Dn); //Precisa corrigir Dw varia no AMR
	  lw.insert(lw.end(), wn);
	}
      }
      else{
	Dp -= Dn;
      }
      wp.set_weight_index((*it)->get_cell_index());
      wp.set_weight_w(Dp);
      lw.insert(lw.end(), wp);
      
      (*it)->set_cell_wpoisson(lw);
      
      //printf("Celula:\n");
      //(*it)->print_cell();
      //printf("pesos:\n");
      //MP = (*it)->get_cell_wpoisson();
      //M->print_weight_list(MP);
      lw.erase(lw.begin(),lw.end());

    }
    
  }
  
  

  
  PetscCall(VecAssemblyBegin(uexa));
  PetscCall(VecAssemblyEnd(uexa));
  PetscCall(VecAssemblyBegin(rhs));
  PetscCall(VecAssemblyEnd(rhs));
  
  /* Create Matrix */
  
  PetscCall(MatCreate(PETSC_COMM_SELF, &A));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, ncell, ncell));
  PetscCall(MatSetFromOptions(A));
  PetscCall(MatSetUp(A));
  
  
  for (i = 0; i < number_of_levels; i++) {
    l = M->get_list_cell_by_level(i);
    
    for (list <cell *>::iterator it = l->begin(); it != l->end(); it++) {
      iglobal = (*it)->get_cell_index();
     
      MP = (*it)->get_cell_wpoisson();
      printf("Celula:\n");
      (*it)->print_cell();
      printf("pesos:\n");
      M->print_weight_list(MP);
      j = 0;
      //printf("%d %d\n",iglobal, j);
      for(list <weight>::iterator wt = MP.begin(); wt != MP.end(); wt++){
	col2[j] = (wt)->get_weight_index();
	w2[j] = (wt)->get_weight_w();
	//printf("%d %d %d\n", iglobal, j, col2[j]);
	j++;
      }
      PetscCall(MatSetValues(A, 1, &iglobal, j, col2, w2, INSERT_VALUES));
    }
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


