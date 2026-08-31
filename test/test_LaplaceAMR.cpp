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
/*
PetscErrorCode ExportToVTK(DM da, Vec u, const char *filename)
{
    PetscViewer viewer;
    PetscErrorCode ierr;

    // 1. Vincula o DMDA ao vetor
    ierr = VecSetDM(u, da); CHKERRQ(ierr);

    // 2. Cria o viewer PETSc
    ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); CHKERRQ(ierr);

    // 3. Define o tipo para VTK (compatível com o VisIt)
    ierr = PetscViewerSetType(viewer, PETSCVIEWERVTK); CHKERRQ(ierr);
    ierr = PetscViewerFileSetMode(viewer, FILE_MODE_WRITE); CHKERRQ(ierr);
    ierr = PetscViewerFileSetName(viewer, filename); CHKERRQ(ierr);

    // 4. Grava os dados da solução
    ierr = VecView(u, viewer); CHKERRQ(ierr);

    // 5. Libera a memória do viewer
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    return 0;
}

PetscErrorCode ExportMeshToVTK(Vec u, int nx, int ny, double lx, double ly, const char *filename)
{
    PetscErrorCode ierr;
    const PetscScalar *u_array;

    // 1. Obtém ponteiro de leitura para o vetor PETSc
    ierr = VecGetArrayRead(u, &u_array); CHKERRQ(ierr);

    // 2. Abre o arquivo ASCII para escrita
    FILE *f = fopen(filename, "w");
    if (!f) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Não foi possível abrir o arquivo VTK.");

    // Cabeçalho padrão do VTK
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Solucao Laplace sem DMDA PETSc\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET STRUCTURED_GRID\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", nx, ny);

    // 3. Escreve as coordenadas dos pontos (x, y, z)
    fprintf(f, "POINTS %d double\n", nx * ny);
    double dx = lx / (nx - 1);
    double dy = ly / (ny - 1);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(f, "%f %f 0.0\n", i * dx, j * dy);
        }
    }

    // 4. Escreve os dados do vetor PETSc associados a cada ponto
    fprintf(f, "\nPOINT_DATA %d\n", nx * ny);
    fprintf(f, "SCALARS solucao double 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < nx * ny; i++) {
        fprintf(f, "%f\n", (double)PetscRealPart(u_array[i]));
    }

    fclose(f);

    // 5. Libera o ponteiro do vetor PETSc
    ierr = VecRestoreArrayRead(u, &u_array); CHKERRQ(ierr);

    return 0;
}
*/
static char help[] = "Solve a Poisson 2D in uniform mesh!";

int main (int argc, char **args){

  int number_of_levels = 3, nxb = 4, nyb;
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
  PetscReal norm, aux, auxp;
  Vec uexa, rhs, r, u;
  Mat A;
  PetscInt i, j, iglobal, N;
  PC pc;

  PetscMPIInt size;
  PetscLogDouble t1, t2;

  /* Geração da malha adaptativa */
  for(int i = 0; i < number_of_levels - 1; i++){
    l = M->get_list_cell_by_level(i);
    dx = fabs(xend - xbegin)/(nxb * pow(2, i));
    dy = fabs(yend - ybegin)/(nyb * pow(2, i));

    it = l->begin();

    while(it != l->end()){
      xd = xbegin + (((*it)->get_cell_x()) * dx);
      yd = ybegin + (((*it)->get_cell_y()) * dy);
      if (xd >= 0.5 && xd < 0.75 && yd >= 0.25 && yd < 0.75){
        if(i == 0){
          it = M->split(*it);
        }
        else
            it++;

      }
      else if( xd >= 0.25 && xd < 0.5 && yd >= 0.5 && yd < 0.75){
	     it = M->split(*it);
      }
      else
	it++;

    }

  }

  /* Enumeração e lista de vizinhas */
  ncell = M->counting_mesh_cell();
  //ncell = M->counting_mesh_cell_inverse();
  N = ncell;
  printf("%d\n", ncell);
  //M->print_mesh();
  //int ncellv = M->neighbours_all_cell();

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
      //printf("%d\n", iglobal);

      // Solucao exata
      aux = (PetscScalar)(sol_u(xd, yd, 0));
      VecSetValues(uexa, 1, &iglobal, &aux, INSERT_VALUES);
      // Lado direito
      aux = dx*dy * (PetscScalar)(f_rhs(xd, yd));
      VecSetValues(rhs, 1, &iglobal, &aux, INSERT_VALUES);

      ppc = M->search(i, j, level);
      cout << i << " " << j << " " << iglobal << endl;

      if (ppc != NULL) {
        lv = M->neighbours(ppc);
        auxp = 4.0;
        idx.push_back(iglobal);

        // Condição de Contorno Dirichlet - variável no centro
        if(i == 0){
           aux = 2.0 * (PetscScalar)(sol_u(xbegin, yd, 0));
           VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
           auxp = auxp + 1.0;
        }
        if(i == nxb * pow(2, kl)-1){
          aux = 2.0 * (PetscScalar)(sol_u(xend, yd, 0));
          VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
          auxp = auxp + 1.0;
        }
        if(j == 0){
          aux = 2.0 * (PetscScalar)(sol_u(xd, ybegin, 0));
          VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
          auxp = auxp + 1.0;
        }
        if(j == nyb * pow(2, kl) - 1){
          aux = 2.0 * (PetscScalar)(sol_u(xd, yend, 0));
          VecSetValues(rhs, 1, &iglobal, &aux, ADD_VALUES);
          auxp = auxp + 1.0;
        }

        wp.push_back(auxp);
        for (list <cell *>::iterator ite = lv->begin(); ite != lv->end(); ite++) {
          idv = (*ite)->get_cell_index();
          idx.push_back(idv);
          wp.push_back(-1.0);
          printf("%d ", idv);
          printf("\n");
        }
        //std::cout << "Tamanho: " << idx.size() << " " << wp.size() << "\n";

        //std::cout << "Vazio?: " << lv->size() << "\n";
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
   
  PetscCall(VecAssemblyBegin(uexa));
  PetscCall(VecAssemblyEnd(uexa));
  PetscCall(VecAssemblyBegin(rhs));
  PetscCall(VecAssemblyEnd(rhs));

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
   //ExportToVTK(da, u, "Resultado_MUE2.vts")
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

   M->create_unstructured_mesh(&sol_u, 0.0);

   return 0;
 }


