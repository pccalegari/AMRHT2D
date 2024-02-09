#include "particle.h"

particle::particle(){
  this->x = this->y = this->vx = this->vy = this->level = this->radius = -1;
}

particle::particle(double x, double y, double radius, double velocity, int level){
  this->x = x;
  this->y = y;
  this->radius = radius;
  //  this->vx = vx;
  //this->vy = vy;
  this->level = level;
}

double particle::get_particle_x(){
  return this->x;
}

double particle::get_particle_y(){
  return this->y;
}

double particle::get_particle_vx(){
  return this->vx;
}

double particle::get_particle_vy(){
  return this->vy;
}

double particle::get_particle_radius(){
  return this->radius;
}

int particle::get_particle_level(){
  return this->level;
}

void particle::set_particle_x(double x){
  this->x = x;
}

void particle::set_particle_y(double y){
  this->y = y;
}

void particle::set_particle_radius(double r){
  this->radius = r;
}

void particle::set_particle_vx(double u){
  this->vx = u;
}

void particle::set_particle_vy(double v){
  this->vy = v;
}

void particle::set_particle_level(int l){
  this->level = l;
}

void particle::print_particle(){
  printf("Coordenada do centro da partícula: (%.4lf, %.4lf)\nRaio: %.4lf\nVelocidade: (%4lf,%.4lf)\nNível da célula que contém o centro da partícula: %d\n",
	 this->x, this->y, this->radius, this->vx, this->vy, this->level);
}
