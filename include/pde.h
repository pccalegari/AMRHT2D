#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

class pde{
  mesh * M;
  double (*bc_north)(double x, double y, double t);
  double (*bc_south)(double x, double y, double t);
  double (*bc_east)(double x, double y, double t);
  double (*bc_west)(double x, double y, double t);

public:
  pde();
  pde(mesh * M, double (*fx1)(double, double, double), double (*fx2)(double, double, double), double (*fy1)(double, double, double), double (*fy2)(double, double, double));
  double dbc_phi(cell * c, double t);
  double nbc_phi(cell * c, double t);
};
