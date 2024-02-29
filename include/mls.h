#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

double delta_s2(double x);
double delta_s1(double x);
double delta_s3(double x);
double delta(double x, int op);
double delta_inter(double x, double y, vector <double> xv, vector <double> yv, vector <double> phi, double d, int op);
double max_var(vector <double> var);
//void pesos_inter_diag(point p, cell c);
