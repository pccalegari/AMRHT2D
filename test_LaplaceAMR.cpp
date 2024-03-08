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

void print_matrix(double **A, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      printf("%f\t", A[i][j]);
    }
    printf("\n");
  }
}

double ** multiplicaAB(double ** A, double **B, int m, int n){
  double ** U;
  U = (double **) malloc(m * sizeof(double *));
  for(int i = 0; i < m; i++)
    U[i] = (double *) malloc(n * sizeof(double));
  double soma;
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      soma = 0.0;
      for(int k = 0; k < n; k++){
	soma += A[i][k]*B[k][j];
      }
      U[i][j] = soma;
    }
  }
  return U;
}

double distance(double px, double py, double qx, double qy){
  double d = sqrt((px - qx)*(px - qx) + (py - qy)*(py - qy));
  return d;
}

double delta_linear(double r){
  double delta = 0.0;
  r = fabs(r);
  if((r >= 0) && (r <= 1))
    delta = 1 - r;
  return delta;
}

double delta_inter(double x, double y, list <cell *> * lv, mesh *M){
  double peso = 0.0, dx, dy, xv, yv, xb, xe, yb, ye, deltax, deltay;
  int nl, i, j, nxb, nyb;
  dominio * D;
  D = M->get_dominio();  
  xb = D->get_xbegin();
  xe = D->get_xend();
  yb = D->get_ybegin();
  ye = D->get_yend();
  nxb = M->get_nxb();
  nyb = M->get_nyb();
  for(list <cell *>::iterator it = lv->begin(); it != lv->end(); it++){
    nl = (*it)->get_cell_level();
    i = (*it)->get_cell_x();
    j = (*it)->get_cell_y();
    dx = fabs(xe - xb)/(nxb * pow(2, nl));
    dy = fabs(ye - yb)/(nyb * pow(2, nl));
    xv = xb + (i + 0.5)*dx;
    yv = yb + (j + 0.5)*dy;
    deltax = delta_linear((xv - x)/dx);
    deltay = delta_linear((yv - y)/dy);
    peso = peso + deltax*deltay*dx*dy;
  }
  
  return peso;
}

double delta_s2(double r){
  double delta = 0.0;
  r = fabs(r);
  if(r >= 0.0 && r < 1.0)
    delta = (1.0/8.0)*(3.0 - 2.0*r + sqrt(1.0 + 4.0*r - 4.0*r*r));
  else if(r >= 1.0 && r <= 2.0)
    delta = (1.0/8.0)*(5.0 - 2.0*r - sqrt(12.0*r - 7.0 - 4.0*r*r));
  return delta;
}

double weight_mls(double x, double y, double xp, double yp){
  double r = sqrt((x - xp)*(x - xp) +(y - yp)*(y - yp));
  double w;
  w = 0;
  if(r <= 0.5){
    w = 2.0/3.0 - 4*r*r + 4*r*r*r;
  }
  if(r > 0.5 && r <= 1.0){
    w = 4.0/3.0 - 4*r + 4*r*r - (4.0/3.0)*r*r*r;
  }
  return w;
}

void double w_delta(double x, double y){
  //localiza cell que contem (x,y)
  //lista de vizinhas da célula
  int nviz = 0;
  // dx e dy da célula que contem (xp,yp)
  
  for(int i = 0; i < nviz; i++){
    //deltax = delta_linear((xv[i] - x)/dx);
    //deltay = delta_linear((yv[i] - y)/dy);
    //wi = deltax*deltay*dx*dy;
    //total += wi;
    
  }
  
  
}

double prod_esc(double * x, double * y, int n){
  double pi = 0.0;
  for(int i = 0; i < n; i++){
    pi += x[i]*y[i];
  }
  return pi;
}
/*
double ** matriz_Q(double v, int m, int k){
  double ** Q;
  double pi = prod_esc(v, v, m-k);
  Q = (double **)malloc((m-k)*sizeof(double *));
  for(int i = 0; i < m-k; i++)
    Q[i] = (double *)malloc((m-k)*sizeof(double));
  for(int i = 0; i < m-k; i++){
    for(int j = 0; j < m-k; j++){
      Q[i][j] = -2.0*v[i]*v[j]/pi;
      if(i == j)
	Q[i][j] += 1.0;
    }
  }
  return Q;
}
*/
double house(double * x, double * v, int n){
  double sigma, beta, mu;
 
  v[0] = 0.0;
  for(int i = 1; i < n; i++)
    v[i] = x[i];
  sigma = prod_esc(v, v, n);
  v[0] = 1.0;
  if (sigma < 1.0e-10)
    beta = 0.0;
  else{
    mu = sqrt(x[0]*x[0] + sigma);
    if (x[0] < 1.0e-10)
      v[0] = x[0] - mu;
    else
      v[0] = -sigma/(x[0] + mu);
    beta = 2.0*v[0]*v[0]/(sigma + v[0]*v[0]);
    for(int i = 0; i < n; i++)
      v[i] = v[i]/v[0];
  }
  return beta;
}

double * peso_mls(double ** Q, double * W, double * y, int m, int n){
  double * mls;
  mls = (double *)malloc(m * sizeof(double));
  double soma;
  
  for(int i = 0; i < m; i++){
    soma = 0.0;
    for(int j = 0; j < n; j++)
      soma += Q[i][j]*W[i]*y[j];
    mls[i] = soma;
  }
  return mls;
}


void solve_RTyPhi(double ** A, double * b, double * y, int n){
  double soma;
  for(int i = 0; i < n; i++){
    soma = 0.0;
    for(int j = 0; j < i; j++)
      soma += A[j][i]*y[j];
    y[i] = (b[i] - soma)/A[i][i];
  }
}

void qr_house(double ** A, int m, int n){
  int tam, i, j, k, l;
  double pi, soma;
  double * x, * v, * beta;
  beta = (double *)malloc(n*sizeof(double));
  for(j = 0; j < n; j++){
    tam = m - j;
    v = (double *)malloc(tam*sizeof(double));
    x = (double *)malloc(tam*sizeof(double));
    for(i = j; i < m; i++)
      x[i] = A[i][j];
    beta[j] = house(x, v, tam);
    pi = prod_esc(v, v, tam);
    for(i = j; i < m; i++){
      for(k = j; k < n; k++){
	soma = 0.0;
	for(l = j; l < m; l++)
	  soma -= beta[j]*(v[k-j]*v[l-j])*A[l][k];
	soma += A[k][k];
	A[i][k] = soma;
      }
    }
    if(j < m){
      for(i = j+1; i < m; i++)
	A[i][j] = v[i-j];
    }
    //printf("%f\n", beta);
    for(i = 0; i < tam; i++)
      printf("%f\t", v[i]);
    printf("\n");
    free(x);
    free(v);
  }
  free(beta);
}

void qr(double ** A, int m, int n, double **Q, double **R){
  for (int j = 0; j < n; j++) {
    // Inicialização
    for (int i = 0; i < m; i++) {
      Q[i][j] = A[i][j];
    }
    for (int k = 0; k < j; k++) {
      // Projeção e subtração
      double dot_product = 0.0;
      for (int i = 0; i < m; i++) {
	dot_product += Q[i][j] * Q[i][k];
      }
      for (int i = 0; i < m; i++) {
	Q[i][j] -= dot_product * Q[i][k];
      }
    }
    // Normalização
    double norm = 0.0;
    for (int i = 0; i < m; i++) {
      norm += Q[i][j] * Q[i][j];
    }
    norm = sqrt(norm);
    for (int i = 0; i < m; i++) {
      Q[i][j] /= norm;
    }
    // Preenchimento da matriz R
    for (int i = 0; i <= j; i++) {
      R[i][j] = 0.0;
      for (int k = 0; k < m; k++) {
	R[i][j] += Q[k][i] * A[k][j];
      }
    }
  }
  
}

void peso(double xp, double yp, list <cell *> *lv, mesh * M, int * index, double * mls){
  double * w, ** A, ** Q, **QQ, ** R, ** RR, ** B, * phi, * ya;
  int m, k, nl, i, j, nxb, nyb, n;
  //n tamanho da base: n = 3 temos {1, x, y}; n=6 temos {1, x, y, x*x, x*y, y*y}
  double dx, dy, xv, yv, xe, ye, xb, yb;
  dominio * D;
  Mat A1;
  k = 0;
  n = 3;
  m = lv->size();
  w = (double *)malloc(m * sizeof(double));
  phi = (double *)malloc(n * sizeof(double));
  ya = (double *)malloc(n * sizeof(double));
  A = (double **) malloc (m * sizeof(double *));
  for(int i = 0; i < m; i++)
    A[i] = (double *) malloc(n * sizeof(double));
  //Q = (double **) malloc (m * sizeof(double *));
  //for(int i = 0; i < m; i++)
  //  Q[i] = (double *) malloc(m * sizeof(double));
  //R = (double **) malloc (m * sizeof(double *));
  //for(int i = 0; i < m; i++)
  //  R[i] = (double *) malloc(n * sizeof(double));
  //QQ = (double **) malloc (m * sizeof(double *));
  //for(int i = 0; i < m; i++)
  //  QQ[i] = (double *) malloc(n * sizeof(double));
  //RR = (double **) malloc (n * sizeof(double *));
  //for(int i = 0; i < n; i++)
  //  RR[i] = (double *) malloc(n * sizeof(double));
  
  D = M->get_dominio();  
  xb = D->get_xbegin();
  xe = D->get_xend();
  yb = D->get_ybegin();
  ye = D->get_yend();
  nxb = M->get_nxb();
  nyb = M->get_nyb();
  list <cell *>::iterator it = lv->begin();
  double h = max(fabs(xe - xb)/(nxb * pow(2, (*it)->get_cell_level())), fabs(ye - yb)/(nyb * pow(2, (*it)->get_cell_level())));
  for(it = lv->begin(); it != lv->end(); it++){
    nl = (*it)->get_cell_level();
    i = (*it)->get_cell_x();
    j = (*it)->get_cell_y();
    
    dy = fabs(ye - yb)/(nyb * pow(2, nl));
    dx = fabs(xe - xb)/(nxb * pow(2, nl));
    dy = fabs(ye - yb)/(nyb * pow(2, nl));
    xv = xb + (i + 0.5)*dx;
    yv = yb + (j + 0.5)*dy;
    double r = sqrt((xv - xp)*(xv - xp) +(yv - yp)*(yv - yp))/h;
    w[k] = delta_s2(r);//1.0;//weight_mls(xp, yp, xv, yv);
    index[k] = (*it)->get_cell_index();
    //matriz WA
    A[k][0] = w[k];
    A[k][1] = w[k]*xv;
    A[k][2] = w[k]*yv;
    k++;
  }
  printf("A:\n");
  print_matrix(A, m, n);
  
  PetscCall(MatCreate(PETSC_COMM_SELF, &A1));
  PetscCall(MatSetSizes(A1, m, n, M, N));
  PetscCall(MatSetFromOptions(A1));
  PetscCall(MatSetUp(A1));
    
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++){
      PetscCall(MatSetValues(A1, m, &i, n, j,(PetscScalar) A[i][j], INSERT_VALUES));
    }
  }
  
  PetscCall(MatAssemblyBegin(A1, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A1, MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(A1, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
  //MatQRFactor(A1, IS columnpermutation, const MatFactorInfo *info);
  //qr(A, m, n, Q, R);
  //B = multiplicaAB(Q, R, m, n);
  //printf("Q:\n");
  //print_matrix(Q, m, m);
  //printf("R:\n");
  //print_matrix(R, m, n);
  //printf("QR:\n");
  //print_matrix(B, m, n);

  //qr_house(A, m, n);

  printf("A = QR:\n");
  print_matrix(A, m, n);
  //for(k = 0; k < m; k++)
  //  for(int l = 0; l < n; l++)
  //    QQ[k][l] = Q[k][l];
  //for(k = 0; k < n; k++)
  //  for(int l = 0; l < n; l++)
  //    RR[l][k] = R[k][l];
  //printf("QQ:\n");
  //print_matrix(QQ, m, n);
  //printf("RR:\n");
  //print_matrix(RR, n, n);
 
  
  //Phi
  phi[0] = 1.0;
  phi[1] = xp;
  phi[2] = yp;
  //printf("%f %f %f\n", phi[0], phi[1], phi[2]);
  //ya = subst_direta(RR, phi, n);
  //printf("Y: %f %f %f\n", ya[0], ya[1], ya[2]);
  //mls = peso_mls(QQ, w, ya, m, n);
  //for(int i = 0; i < m; i++)
  //  printf("Pesos: %f\t", mls[i]);
  //printf("\n");
  solve_RTyPhi(A, phi, ya, n);
  printf("Y:\n");
  for(int i = 0; i < n; i++)
    printf("%f\t", ya[i]);
  printf("\n");
   for(int i = 0; i < m; i++){
    //free(Q[i]);
    //free(R[i]);
    free(A[i]);
  }
  //free(Q);
  //free(R);
  free(A);
  //for(int i = 0; i < m; i++)
  //  free(QQ[i]);
  //for(int i = 0; i < n; i++)
  //  free(RR[i]);

  
  //free(QQ);
  //free(RR);

  free(w);
  free(phi);
  free(ya);
  
}
/*
weight busca(int x, list <weight> lw){
  list <weight>::iterator k = lw.begin();
  weight ww;
  while(k != lw.end() && k->get_weight_index() != x)
    k++;
  if(k != lw.end())
    ww = NULL;
  else
    ww = k;
  return ww;
}
*/
/*
void insere_peso_mls(int * index, double * mls, int m, list <weight> * lw) {
  weight ww;
  int tam = lw->size();
  
  for(int i = 0; i < m; i++){
    ww.set_weight_index(index[k]);
  ww.set_weight_w(Dw); 
  lw.insert(lw.end(), ww);
}

*/

//extern PetscErrorCode ComputeMatrix(KSP, Mat, Mat, void *);
//extern PetscErrorCode ComputeRHS(KSP, Vec, void *);
//extern PetscErrorCode ComputeInitialGuess(KSP, Vec, void *);

static char help[] = "Solve a Poisson 2D!";

int main (int argc, char **args){

  int number_of_levels = 2, nxb = 4, nyb = 4;
  int ncell, m, nivel, nl;
  dominio *D;
  double xbegin, ybegin, xend, yend, xe, xw, yn, ys;
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
  list <cell *> * lvw, * lve, * lvn, * lvs, lvv;
  list <cell *>::iterator it, kt;
  double w2[5];
  double * mls;
  int * id;
  
  KSP ksp;
  PetscReal norm, aux;
  Vec uexa, u, rhs, r;
  Mat A;
  PetscInt i, j, k, iglobal, col2[5];
  PC pc;
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
      if (xd >= 0.5 && xd < 0.75 && yd >= 0.25 && yd < 0.75){
	it = M->split(*it);
      }
      else if( xd >= 0.25 && xd < 0.5 && yd >= 0.5 && yd < 0.75){
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
      (*it)->set_cell_lv(M->neighbours(*it));
      //printf("Celula:\n");
      //(*it)->print_cell();
      //printf("Vizinhas:\n");
      //M->print_cell_list((*it)->get_cell_lv());
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
    aux = - Ds * 2.0 * (PetscScalar)(sol_u(xd, ybegin));
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
    aux = - Dn * 2.0 * (PetscScalar)(sol_u(xd, yend));
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
    aux = - De * 2.0 * (PetscScalar)(sol_u(xend, yd));
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
    aux = - Dw * 2.0 * (PetscScalar)(sol_u(xbegin, yd));
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
      printf("Celula:\n");
      (*it)->print_cell();
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      xd = xbegin + ((i + 0.5) * dx);
      yd = ybegin + ((j + 0.5) * dy);
      xe = xd + dx;
      xw = xd - dx;
      yn = yd + dy;
      ys = yd - dy;
      De = -Dif(xd + dx*0.5, yd)/dx2;
      Dw = -Dif(xd - dx*0.5, yd)/dx2;
      Ds = -Dif(xd, yd - dy*0.5)/dy2;
      Dn = -Dif(xd, yd + dy*0.5)/dy2;
      lvw = (*it)->get_cell_lbw();
      lvs = (*it)->get_cell_lbs();
      lve = (*it)->get_cell_lbe();
      lvn = (*it)->get_cell_lbn();
      lvv.assign(lvw->begin(), lvw->end());
      //lv->insert(lv->end(), lvw->begin(), lvw->end());
      lvv.insert(lvv.end(), lvs->begin(), lvs->end());
      lvv.insert(lvv.end(), lve->begin(), lve->end());
      lvv.insert(lvv.end(), lvn->begin(), lvn->end());
      lvv.insert(lvv.end(),(*it));
     
      //M->print_cell_list(&lvv);
      Dp = -(De + Dw + Ds + Dn);
      //printf("Tamanhos %ld %ld %ld %ld\n", lvw->size(), lvs->size(), lve->size(), lvn->size());
      if((*it)->get_cell_nvw() > 1){
	//printf("Vizinhas finas: W\n");
	m = ((*it)->get_cell_lv())->size();
	mls = (double *)malloc(m * sizeof(double));
	id = (int *)malloc(m * sizeof(int));
	peso(xw, yd, (*it)->get_cell_lv(), M, id, mls);
	for(int a = 0; a < m; a++)
	  printf("%d %f\n", id[a], mls[a]);
	//Vizinhas para o MMQ
	//int k = 0;
	//for(list <cell *>::iterator wt = lv.begin(); wt != lv.end(); wt++){
	//  weight idx = busca(id[k], lw);
	//  printf("%d %d %f\n", (*wt)->get_cell_index(), idx.get_weight_index(), idx.get_weight_w());
	  //if(idx < 0){
	    //Atualiza pesos
	    //ww.set_weight_index(id[k]);
	    //   ww.set_weight_w(mls[k]); 
	    //lw.insert(lw.end(), ww);
	    //k++;
	  //}
	  //else{
	    
	  //}
	//}
      }
      else if((*it)->get_cell_nvw() < 1){
	//printf("Vizinha fronteira: W\n");
	Dp -= Dw;
      }
      else{
	list <cell *>::iterator wt = lvw->begin();
	nl = (*wt)->get_cell_level();
	if(nivel == nl){
	  //printf("Vizinha mesmo nível: W\n");
	  ww.set_weight_index((*wt)->get_cell_index());
	  ww.set_weight_w(Dw);
	  lw.insert(lw.end(), ww);	  
	}
	else{
	  //(*it)->print_cell();
	  //printf("Vizinha mais grossa: W\n");
	}
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
      lvv.erase(lvv.begin(),lvv.end());

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
      //printf("Celula:\n");
      //(*it)->print_cell();
      //printf("pesos:\n");
      //M->print_weight_list(MP);
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


