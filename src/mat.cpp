#include <cmath>
#include <cassert>
#include "mesh.h"
#include "mat.h"

mat::mat(){
  A = NULL;
  n = -1;
  IA = NULL;
  JA = NULL;
}

/*void mesh::rhs_dirichlet_boundary_conditions (double (*f) (double x, double y), vector<double> * fvalue) {
  double xbegin, xend, ybegin, yend, xmiddle, ymiddle;
  double dx, dy, dx2, dy2;
  int i, j;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();
  for (int ll = 0; ll < number_of_levels; ll++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, ll));
    dx2 = dx*dx;
    dy = fabs(yend - ybegin) / (nyb * pow(2, ll));
    dy2 = dy*dy;
    for (list <cell *>::iterator it = l->at(ll)->begin(); it != l->at(ll)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      xmiddle = xbegin + (i + 0.5) * dx;
      ymiddle = ybegin + (j + 0.5) * dy;

      cout << ll << " " << i << " " << j << " " << xmiddle << " " << ymiddle << " " << xbegin << " " << xend << " " << dx << endl;
      
      if(j == 0)
	fvalue->push_back((*f) (xmiddle, ymiddle) + 2.0 * (*f) (xmiddle, ybegin)/dy2);
      
      else if (j == pow(2,ll)*nyb - 1)
	fvalue->push_back((*f) (xmiddle, ymiddle) + 2.0 * (*f) (xmiddle, yend)/dy2);
	
      else if(i == 0)
	fvalue->push_back((*f) (xmiddle, ymiddle) + 2.0 * (*f) (xbegin, ymiddle)/dx2);
      
      else if (i == pow(2,ll)*nxb - 1)
	fvalue->push_back((*f) (xmiddle, ymiddle) + 2.0 * (*f) (xend, ymiddle)/dx2);
      else
	fvalue->push_back((*f) (xmiddle, ymiddle));
    }
  }
  } */

