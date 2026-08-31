#include <cstdio>
#include <cstdlib>

using namespace std;

class particle {
 private:
  double x;
  double y;
  double vx;
  double vy;
  double radius;
  int level;
  
 public:
  particle();
  particle(double x, double y, double r, double velocity, int level);
  double get_particle_x();
  double get_particle_y();
  double get_particle_vx();
  double get_particle_vy();
  double get_particle_radius();
  int get_particle_level();
  void set_particle_x(double x);
  void set_particle_y(double y);
  void set_particle_radius(double r);
  void set_particle_vx(double u);
  void set_particle_vy(double v);
  void set_particle_level(int level);
  void print_particle();
};
