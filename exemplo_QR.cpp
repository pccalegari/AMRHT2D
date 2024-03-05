#include <cmath>
#include <stdio.h>
#include <iostream>

void print_matrix(double **A, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      printf("%f\t", A[i][j]);
    }
    printf("\n");
  }
}

double ** transpose(double **A, int m, int n){
  double **B;
  B = (double **) malloc (n * sizeof(double *));
  for(int i = 0; i < n; i++)
    B[i] = (double *) malloc(m * sizeof(double));
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      B[j][i] = A[i][j];
    }
  }
  return B;
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

double ** prod_v_vt(double * v, double beta, int m, int n){
  double aux = 0;
  double ** U;
  U = (double **) malloc(m * sizeof(double *));
  for(int i = 0; i < m; i++)
    U[i] = (double *) malloc(n * sizeof(double));
  for(int i = 0; i < m; i++){
    for(int j = 0; j < m; j++){
      if(i == j)
	aux = 1;
      U[i][j] = aux - beta*v[i]*v[j];
    }
  }
  return U;
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
      for(int k = 0; k < m; k++){
	soma += A[i][k]*B[k][j];
      }
      U[i][j] = soma;
    }
  }
  return U;
}

double norma(double * v, int n){
  double aux = 0.0, norma;
  for(int i = 0; i < n; i++)
    aux += v[i]*v[i];
  norma = sqrt(aux);
  return norma;
}

double prod_escalar(double * v, double *u, int n){
  double aux = 0.0;
  for(int i = 0; i < n; i++)
    aux += v[i]*u[i];
  return aux;
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


void house(double * v, double beta, double * x, int n){
  double sigma, mu;
  sigma = 0;
  for(int i = 1; i < n; i++){
    sigma += x[i]*x[i];
  }
  v[0] = 1;
  for(int i = 1; i < n; i++){
    v[i] = x[i];
  }
  if(sigma == 0)
    beta = 0;
  else{
    mu = sqrt(x[0]*x[0] + sigma);
    if(x[0] <= 0){
      v[0] = x[0] - mu;
    }
    else{
      v[0] = -sigma/(x[0] + mu);
    }
    beta = (2*v[0]*v[0])/(sigma + v[0]*v[0]);
    for(int i = 0; i < n; i++)
      v[i] = v[i]/v[0];
  }
    
}

int main(int argc, char **argv) {
  int m = 6, n = 3;
  double x = 0.0;
  double **A, **B, **Q, **R, **C;
  A = (double **) malloc (m * sizeof(double *));
  for(int i = 0; i < m; i++)
    A[i] = (double *) malloc(n * sizeof(double));
  C = (double **) malloc (m * sizeof(double *));
  for(int i = 0; i < m; i++)
    C[i] = (double *) malloc(n * sizeof(double));
  Q = (double **) malloc (m * sizeof(double *));
  for(int i = 0; i < m; i++)
    Q[i] = (double *) malloc(m * sizeof(double));
  R = (double **) malloc (m * sizeof(double *));
  for(int i = 0; i < m; i++)
    R[i] = (double *) malloc(n * sizeof(double));
  for(int i = 0; i < m; i++){
    A[i][0] = 1.0;
    A[i][1] = x;
    A[i][2] = x*x;
    x += 0.2;
  }
  //A[0].assign(l0,l0+3);
  print_matrix(A, m, n);
  B = transpose(A, m, n);
  printf("\n");
  qr(A, m, n, Q, R);
  print_matrix(R, m, n);
  printf("\n");
  print_matrix(Q, m, m);

  C = multiplicaAB(Q, R, m, n);
  printf("\n");
  print_matrix(C, m, n);

  
  for(int i = 0; i < m; i++){
    free(Q[i]);
    free(R[i]);
    free(A[i]);
    free(C[i]);
  }
  free(Q);
  free(R);
  free(A);
  free(C);
  return 0;
}
