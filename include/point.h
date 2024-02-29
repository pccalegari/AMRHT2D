#include <cstdio>
#include <cstdlib>

using namespace std;

class point {
 private:
  double x, y;
    
 public:
  point ();
  point (double x, double y);
  double get_point_x();
  double get_point_y();
  void set_point_x(double xp);
  void set_point_y(double py);
  double dist(point p, point q);
  //double dist(double px, double py, double qx, double qy);
};
