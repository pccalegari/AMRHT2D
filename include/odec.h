#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>

using namespace std;

class ode{
  double tini;          // tempo inicial
  double sini;          // Condicao Inicial
  double tend;
  double (*sfn)(double t, double phi);  // Funcao fonte

public:
  ode();
  ode(double t0, double phi0, double tend, double (*f)(double, double));
  double rk1(double dt);
  double rk2M(double dt);
  double rk2A(double dt);
  double rk4(double dt);
};
