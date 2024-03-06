#include <cmath>
#include <cassert>
#include "mesh.h"
#include <vector>

mesh::mesh(){
  l = NULL;
  number_of_levels = -1;
  max_dimension_by_level = NULL;
  nxb = -1;
  nyb = -1;
}

mesh::mesh(dominio * D, int n_levels, int nx, int ny){
  this->number_of_levels = n_levels;
  this->max_dimension_by_level = new vector<int>;
  
  for(int i = 0; i < n_levels; i++){
    this->max_dimension_by_level->push_back((nx * pow(2,i)) * (ny * pow(2,i)));
  }
  this->D = D;
  //o tamanho da tabela hash vai depender do máximo número de células do nível mais fino, isto é:
  //(1/100)*max_dim_by_level[n_level - 1] => (x * 100) % de max_dim_by_level[n_level - 1]
  double x = 0.1;
  H = new hash_table((int)(x * max_dimension_by_level->at(n_levels - 1)), n_levels, nx * pow(2, n_levels));
  //H = new hash_table(1000000);
  l = new vector<list<cell *> * >;
  for (int i = 0; i < n_levels; i++){
    list<cell *> * l_by_level= new list <cell *>;
    l->push_back(l_by_level); 
  }
  //number_of_levels = n_levels;
  nxb = nx;
  nyb = ny;
  //this->max_dimension_by_level = max_dimension_by_level;
  cell * c;
  //int index = 0;
  //double iniciau = 0.0, iniciav = 0.0, iniciaphi = 0.0, iniciaphi0 = 0.0;

  for (int y = 0; y < nyb; y++)
    for (int x = 0; x < nxb; x++){
      c = new cell(x, y, 0);
      this->insert(c);
    }
}

void mesh::insert(cell * c){
  c->set_cell_pointer_to_list((l->at(c->get_cell_level()))->insert(l->at(c->get_cell_level())->begin(), c));
  H->insert(c); 
}

list <cell *>::iterator mesh::remove(cell * c) {
  list <cell *>::iterator it;
  it = l->at(c->get_cell_level())->erase(c->get_cell_pointer_to_list());
  H->remove(c);
  delete c;
  return it;
}

int mesh::get_nxb(){
  return nxb;
}

int mesh::get_nyb(){
  return nyb;
}

cell * mesh::search(int x, int y, int level) {
  return H->search(x, y, level);
}

cell * mesh::search_particle_cell(particle * p, int level) {
  double xb, xe, yb, ye, dx, dy;
  int ip, jp;
  xb = get_dominio()->get_xbegin();
  yb = get_dominio()->get_ybegin();
  xe = get_dominio()->get_xend();
  ye = get_dominio()->get_yend();
  dx = fabs(xe - xb)/(nxb*pow(2,level));
  dy = fabs(ye - yb)/(nyb*pow(2,level));
  ip = floor((p->get_particle_x() - xb)/dx);
  jp = floor((p->get_particle_y() - yb)/dy);
  return H->search(ip, jp, level);
}

void mesh::print_mesh() {
  for (int i = 0; i < number_of_levels; i++)
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++)
      cout << (*it)->get_cell_x() << " " << (*it)->get_cell_y() << " " << (*it)->get_cell_level() << " " << (*it)->get_cell_index() << endl;
  
}

hash_table * mesh::get_hash_table() {
  return H;
}


list <cell *>::iterator mesh::split(cell * c) {
  list <cell *>::iterator it;
  cell ** V;
  V = c->split();
  it = remove(c);
  insert(V[0]);
  insert(V[1]);
  insert(V[2]);
  insert(V[3]);
  free (V);
  return it;
}

cell ** mesh::split_return_new_cells(cell * c) {
  list <cell *>::iterator it;
  cell ** V;
  V = c->split();
  it = remove(c);
  insert(V[0]);
  insert(V[1]);
  insert(V[2]);
  insert(V[3]);
  return V;
}

void mesh::create_unstructured_mesh(double (* f) (double x, double y, double t), double tempo) {
  //4 pontos para cada celula (haverá pontos para duas ou mais células, com isso, o fator de carga de internal_HT será menor que 1
  //int x, y;
  //int count = 0;
  double xd, yd, dx, dy;
  vector <double> ptx;
  vector <double> pty;
  vector <int> verticeslist;
  //int number_node = 1;
  double_cell * cs [4];
  double_cell * c, * cmiddle;
  //por enquanto as variáveis abaixo (a, b, delta_0) descrevem o domínio. (a,b) é o ponto onde começa o domínio e delta_0 é o fator de discretização no nível mais grosso (nível 0)
  double xbegin, ybegin, xend, yend;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();

  /****************SHAPE OR ZONES OR CELL ARE THE SAME*************************/
  int ncells = H->get_number_cell();
  int * celltypes = (int *) malloc (sizeof(int) * ncells);
  int * cellsizes = (int *) malloc (sizeof(int) * ncells);
  int * cellcnt = (int *) malloc (sizeof(int) * ncells); // I still confuse about the mean of this variable
  double * cellmiddles = (double *) malloc (sizeof(double) * ncells);
  //double * rhs = (double *) malloc (sizeof(double) * ncells);
  /****************************************************************************/
  
  double_hash_table * internal_HT = new double_hash_table(4 * H->get_number_cell());

  int nzones = 0;
  for (int i = 0; i < number_of_levels; i++) {
    //printf ("Nível: %d\n", i);
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      //cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
      xd = xbegin + ((*it)->get_cell_x() * dx);//(delta_0 / pow(2, i)));
      yd = ybegin + ((*it)->get_cell_y() * dy);//(delta_0 / pow(2, i)));
      
     //clock or counter-clock
      cs[0] = new double_cell (xd, yd, dx, dy);
      cs[1] = new double_cell (xd, yd + dy, dx, dy);
      cs[2] = new double_cell (xd + dx, yd + dy, dx, dy);
      cs[3] = new double_cell (xd + dx, yd, dx, dy); 
      cmiddle = new double_cell (xd + (dx / 2.), yd + (dy / 2.), 0., 0.);
      //count++;//conta o número de células

      for (int j = 0; j < 4; j++) {
	c = internal_HT->search(cs[j]->get_double_cell_x(), cs[j]->get_double_cell_y());
	if (c == NULL) {
	  internal_HT->insert(cs[j]);
	  ptx.push_back(cs[j]->get_double_cell_x());
	  pty.push_back(cs[j]->get_double_cell_y());
	  cs[j]->set_double_cell_pointer_to_vector(ptx.size() - 1);
	  c = cs[j];
	}
	else {
	  delete cs[j];
	}
	verticeslist.push_back(c->get_double_cell_pointer_to_vector());
      }
      //add a ZONE related to the four dots cs[0], cs[1], cs[2] and cs[3]
      celltypes[nzones] = DB_ZONETYPE_QUAD;
      cellsizes[nzones] = 4;
      cellcnt[nzones] = 1;
      cellmiddles[nzones] = (*f)(cmiddle->get_double_cell_x(), cmiddle->get_double_cell_y(), tempo);
      //rhs[nzones] = (*df)(cmiddle->get_double_cell_x(), cmiddle->get_double_cell_y());
      delete cmiddle;
      nzones++;
    }
  }

  assert (nzones == ncells);

  /********************silo code*******************************/
  DBfile *dbfile = NULL;
  /* Open the Silo file */
  char palavra [100];
  for (int i = 0; i < 100; i++)
    palavra[0] = '\0';
  sprintf(palavra,"/home/priscila/Documentos/Pesquisa/AMRHT2D-vizinhas-pesos/AMRHT2D-master/data/basic%4.3lf.silo", 1 + tempo);
  
  dbfile = DBCreate(palavra, DB_CLOBBER, DB_LOCAL,"Comment about the data", DB_HDF5);
  if(dbfile == NULL)
    {
      fprintf(stderr, "Could not create Silo file!\n");
      fprintf(stderr, "Change 'palavra' in src/mesh.cpp!\n");  
      exit(0);
    }
  
  /* Add other Silo calls here. */
  int * verticeslistC = &(verticeslist[0]);
  //double * ptxC = &(ptx[0]);
  //double * ptyC = &(pty[0]);
  double * coords [2];
  coords[0] = (double *) malloc (sizeof(double) * ptx.size());
  coords[1] = (double *) malloc (sizeof(double) * ptx.size());
  for (unsigned int i = 0; i < ptx.size(); i++){
    coords[0][i] = ptx[i];
    coords[1][i] = pty[i];
  }
  
  int lverticeslistC = sizeof(int) * verticeslist.size() / sizeof(int);
  
  int nnodes = ptx.size();
  
  int ndims = 2;

  DBPutZonelist2 (dbfile, "zonelist", nzones, ndims, verticeslistC, lverticeslistC, 0, 0, 0, celltypes, cellsizes, cellcnt, ncells, NULL);
  DBPutUcdmesh (dbfile, "meshBLA", ndims, NULL, coords, nnodes, nzones, "zonelist", NULL, DB_DOUBLE, NULL);
  assert(DBPutUcdvar1 (dbfile, "cellmiddles", "meshbBLA", cellmiddles, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);
  assert(DBPutUcdvar1 (dbfile, "f", "meshBLA", cellmiddles, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);

  //assert(DBPutUcdvar1 (dbfile, "rhs", "meshBLA", cellmiddles, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);
      
  /* Close the Silo file. */
  DBClose(dbfile);  
  /************************************************************/

  /********************check repetitions***********************/
  /*vector<double>::iterator itx = ptx.begin();
  vector<double>::iterator ity = pty.begin();
  while (itx != ptx.end()){
    vector<double>::iterator itx2 = ptx.begin();
    vector<double>::iterator ity2 = pty.begin();
    int count = 0;
    while (itx2 != ptx.end()) {
      if ((*itx) == (*itx2) && (*ity) == (*ity2))
	count++;
      itx2++;
      ity2++;
    }
    assert(count == 1);
    itx++;
    ity++;
    }*/
  /************************************************************/
  //PRINT COORDENATES
  /*itx = ptx.begin();
  ity = pty.begin();
  int i = 1;
  while(itx != ptx.end()) { 
    cout << i << ": " << "(" << (*itx) << ", " << (*ity) << ")" << endl;
    itx++;
    ity++;
    i++;
  }
  cout << endl;*/

  //internal_HT->print_information();

  //TODO
  internal_HT->clean();
  
}

void mesh::print_silo(int ct, list <particle *> * P) {
  //4 pontos para cada celula (haverá pontos para duas ou mais células, com isso, o fator de carga de internal_HT será menor que 1
  //int x, y;
  //int count = 0;
  double xd, yd, dx, dy;
  vector <double> ptx;
  vector <double> pty;
  vector <int> vert;
  //int number_node = 1;
  //por enquanto as variáveis abaixo (a, b, delta_0) descrevem o domínio. (a,b) é o ponto onde começa o domínio e delta_0 é o fator de discretização no nível mais grosso (nível 0)
  double xbegin, ybegin, xend, yend;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();
  
  /****************SHAPE OR ZONES OR CELL ARE THE SAME*************************/
  int ncells = H->get_number_cell();
  int * verticeslist = (int *) malloc (sizeof(int) * (4*ncells));
  int * celltypes = (int *) malloc (sizeof(int) * ncells);
  int * cellsizes = (int *) malloc (sizeof(int) * ncells);
  int * cellcnt = (int *) malloc (sizeof(int) * ncells); // I still confuse about the mean of this variable
  double * u = (double *) malloc (sizeof(double) * ncells);
  double * v = (double *) malloc (sizeof(double) * ncells);
  double * phi = (double *) malloc (sizeof(double) * ncells); 
  /****************************************************************************/
  int i, j;
  int nzones = 0;
  int cnodes = 0;
  for (int k = 0; k < number_of_levels; k++) {
    //printf ("Nível: %d\n", i);
    dx = fabs(xend - xbegin) / (nxb * pow(2, k));
    dy = fabs(yend - ybegin) / (nyb * pow(2, k));
    for (list <cell *>::reverse_iterator it = l->at(k)->rbegin(); it != l->at(k)->rend(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      xd = xbegin + (i * dx);
      yd = ybegin + (j * dy);
      (*it)->set_cell_index(nzones);
      //Aqui a inclusao das coordenadas no sentido horario
      ptx.push_back(xd);
      pty.push_back(yd);
      //vert.push_back(cnodes);
      verticeslist[cnodes] = cnodes;
      cnodes++;
      ptx.push_back(xd);
      pty.push_back(yd + dy);
      verticeslist[cnodes] = cnodes;
      ptx.push_back(xd + dx);
      pty.push_back(yd + dy);
      cnodes++;
      verticeslist[cnodes] = cnodes;
      ptx.push_back(xd + dx);
      pty.push_back(yd);
      cnodes++;
      verticeslist[cnodes] = cnodes;
      /*if(i == pow(2,k)*nxb - 1){
	ptx.push_back(xd + dx);
	pty.push_back(yd);
	vert.push_back(cnodes);
	cnodes++;
      }
      if(j == pow(2,k)*nyb - 1){
	ptx.push_back(xd);
	pty.push_back(yd + dy);
	vert.push_back(cnodes);
	cnodes++;
      }
      if((j == pow(2,k)*nyb - 1) && (i == pow(2,k)*nxb - 1)){
	ptx.push_back(xd + dx);
	pty.push_back(yd + dy);
	vert.push_back(cnodes);
	cnodes++;
	}*/
      //add a ZONE related to the four dots cs[0], cs[1], cs[2] and cs[3]
      celltypes[nzones] = DB_ZONETYPE_QUAD;
      cellsizes[nzones] = 4;
      cellcnt[nzones] = 1;
      u[nzones] = (*it)->get_cell_velu();
      v[nzones] = (*it)->get_cell_velv();
      phi[nzones] = (*it)->get_cell_phi();
      nzones++;
      cnodes++;
    }
  }

  assert (nzones == ncells);

  /********************silo code*******************************/
  DBfile *dbfile = NULL;
  /* Open the Silo file */
  char palavra [100];
  for (int i = 0; i < 100; i++)
    palavra[0] = '\0';
  sprintf(palavra,"/home/priscila/Pesquisa/AMRHT2D-master/data/basic%d.silo", ct);
  
  dbfile = DBCreate(palavra, DB_CLOBBER, DB_LOCAL,"Comment about the data", DB_HDF5);
  if(dbfile == NULL)
    {
      fprintf(stderr, "Could not create Silo file!\n");
      exit(0);
    }
  
  int ndims = 2;
  double * coords [2];
  coords[0] = (double *) malloc (sizeof(double) * ptx.size());
  coords[1] = (double *) malloc (sizeof(double) * ptx.size());
  
  for (unsigned int i = 0; i < ptx.size(); i++){
    coords[0][i] = ptx[i];
    coords[1][i] = pty[i];
  }
  
  int nnodes = ptx.size();

  DBPutZonelist2 (dbfile, "zonelist", nzones, ndims, verticeslist, cnodes, 0, 0, 0, celltypes, cellsizes, cellcnt, ncells, NULL);

  DBPutUcdmesh (dbfile, "meshBLA", ndims, NULL, coords, nnodes, nzones, "zonelist", NULL, DB_DOUBLE, NULL);

  assert(DBPutUcdvar1 (dbfile, "velu", "meshBLA", u, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);
  
  assert(DBPutUcdvar1 (dbfile, "velv", "meshBLA" , v, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);

  assert(DBPutUcdvar1 (dbfile, "phi", "meshBLA", phi, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);

  
  //silo particle
  double *x , *y, *up, *vp, *rp;
  int NPTS = P->size();
  x = (double *) malloc(sizeof(double)*NPTS);
  y = (double *) malloc(sizeof(double)*NPTS);
  up = (double *) malloc(sizeof(double)*NPTS);
  vp = (double *) malloc(sizeof(double)*NPTS);
  rp = (double *) malloc(sizeof(double)*NPTS);
  
  double * coords_points [2];
  coords_points [0] = x;
  coords_points [1] = y;
  
  i = 0;
  for (list<particle *>::iterator it = P->begin(); it != P->end(); it++){
    x[i] = (*it)->get_particle_x();
    y[i] = (*it)->get_particle_y();
    up[i] = (*it)->get_particle_vx();
    vp[i] = (*it)->get_particle_vy();
    rp[i] = (*it)->get_particle_radius();
    i++;
  }
  
  DBPutPointmesh(dbfile, "pointmesh", ndims, coords_points, NPTS, DB_DOUBLE, NULL);

  //DBPutPointvar(dbfile, "veup", "pointmesh", NPTS, up, NPTS, DB_DOUBLE, NULL);

  /* Close the Silo file. */
  DBClose(dbfile);  ;
  
}

void mesh::initialize_var(double (* u) (double x, double y, double t), double (* v) (double x, double y, double t), double (*phi) (double x, double y, double t), double tempo, double tempo0){
  double xbegin, ybegin, xend, yend, dx, dy, xd, yd;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();

  for (int i = 0; i < number_of_levels; i++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      //cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
      xd = xbegin + (0.5 + (*it)->get_cell_x()) * dx;//(delta_0 / pow(2, i)));
      yd = ybegin + (0.5 + (*it)->get_cell_y()) * dy;//(delta_0 / pow(2, i)));
      (*it)->set_cell_velu((*u)(xd, yd, tempo));
      (*it)->set_cell_velv((*v)(xd, yd, tempo));
      (*it)->set_cell_phi((*phi)(xd, yd, tempo));
      (*it)->set_cell_phi0((*phi)(xd, yd, tempo0));
      //(*it)->set_cell_with_particle(0);
    }
  }

}

void mesh::calculation_exact_function (double (* f) (double_cell * c)){
  double xbegin, ybegin, xend, yend, dx, dy, xmiddle, ymiddle;
  double_cell * c;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();
  for (int i = 0; i < number_of_levels; i++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, i));
    dy = fabs(yend - ybegin) / (nyb * pow(2, i));
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      xmiddle = xbegin + (((*it)->get_cell_x() + 0.5) * dx);//(delta_0 / pow(2, i)));
      ymiddle = ybegin + (((*it)->get_cell_y() + 0.5) * dy);//(delta_0 / pow(2, i)));
      
      c = new double_cell (xmiddle, ymiddle, dx, dy);
      
      printf ("%.7f ", (*f)(c));
      
      delete c;
    }
    printf ("\n");
  }
}

list <cell *> * mesh::get_list_cell_by_level (int level){
  return l->at(level);
}

dominio * mesh::get_dominio(){
  return D;
}

list <cell *> * mesh::left_finest_or_same_level_neighbours(cell * c, int level_k) {
  list <cell*> * l = new list <cell *>;
  //O valor em level_k DEVE SER maior ou igual que c->get_cell_level();
  int k = level_k - c->get_cell_level();
  assert (k >= 0);
  int number_neighbours_by_x = (int) pow(2, k);
  int number_neighbours_by_y = (int) pow(2, k);
  int _2_to_k = number_neighbours_by_x;
  int i, j, x, y;
  x = c->get_cell_x();
  y = c->get_cell_y();
  cell * neighbour_c;
  for (j = 0; j < number_neighbours_by_y; j++) {
    for (i = 0; i < number_neighbours_by_x; i++) {
      //(x - 1) determina os vizinhos à esquerda de c
      neighbour_c = search (_2_to_k * (x - 1) + i, _2_to_k * y + j, level_k);
      assert(neighbour_c != NULL);
      l->push_back(neighbour_c);
    }
  }
  return l;
}

int mesh::left_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  int n = 0;
  if (x != 0) {
    int l;
    if (x % 2 == 0){
      for (l = level - 1; l >= 0; l--){
	xl = (int) (x / pow (2, level - l));
	yl = (int) (y / pow (2, level - l));
	xl = xl - 1;
	//printf ("(%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, l);
	if (found != NULL) {
	  nb->insert (nb->end(), found);
	  n++;
	  //cout << "FOUND in coarser" << endl;
	  //printf ("L.: (%d, %d): %d)\n", xl, yl, l);
	  break;
	}
      }
    }
    if (x % 2 == 1 || (x % 2 == 0 && l == -1)) {//não encontrou uma célula à esquerda mais grossa
      int xx;
      for (l = level; l < number_of_levels; l++) {
	xx = (int) (pow (2, l - level) * x);
	for (xl = xx - 1; xl < xx; xl++){
	  for (int j = 0; j <= (int) pow (2, l - level) - 1; j++) {
	    yl = (int) (pow (2, l - level) * y) + j;
	    //printf ("(%d, %d): %d)\n", xl, yl, i);
	    cell * found = H->search (xl, yl, l);
	    if (found != NULL) {
	      nb->insert (nb->end(), found);
	      n++;
	      ///cout << "FOUND in finer" << endl;
	      //printf ("L.: (%d, %d): %d)\n", xl, yl, l);
	    }
	  }
	} 
      }
    }
  }
  return n;
  //else{
    
  //}
}

int mesh::right_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  int n = 0;
  if (x != nxb*pow(2,level)) {
    int i;
    if (x % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	xl = xl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->end(), found);
	  n++;
	  //cout << "FOUND in coarser" << endl;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (x % 2 == 0 || (x % 2 == 1 && i == -1)) {//não encontrou uma célula à esquerda mais grossa 
      for (i = level; i < number_of_levels; i++) {
	xl = (int) (pow (2, i - level) * x);
	xl = xl + pow(2, i - level); //procuramos pelos vizinhos à direita de c e a referência de c é o canto inferior esquerdo. Temos que "pular" as células nível atual (+ fino) dentro do próprio c 
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  yl = (int) (pow (2, i - level) * y) + j;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->end(), found);
	    n++;
	    //cout << "FOUND in finer" << endl;
	    //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
  return n;
}

int mesh::down_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  int n = 0;
  if (y != 0) {
    int i;
    if (y % 2 == 0){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl - 1;
	//printf ("(%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL){
	  nb->insert (nb->end(), found);
	  n++;
	  //found = H->search(xl-1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("DL.: (%d, %d): %d)\n", xl-1, yl, i);
	  //}
	  //found = H->search(xl+1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("DR.: (%d, %d): %d)\n", xl+1, yl, i);
	  //}
	  //cout << "FOUND in coarser" << endl;
	  //printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (y % 2 == 1 || (y % 2 == 0 && i == -1)) {//não encontrou uma célula abaixo mais grossa 
      for (i = level; i < number_of_levels; i++) {
	yl = (int) (pow (2, i - level) * y);
	yl = yl - 1;
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  xl = (int) (pow (2, i - level) * x) + j;
	  //printf ("(%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->end(), found);
	    n++;
	    //cout << "FOUND in finer" << endl;
	    //printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	    //found = H->search (xl-1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf ("DL.: (%d, %d): %d)\n", xl-1, yl, i);
	    //}
	    //found = H->search (xl+1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf ("DR.: (%d, %d): %d)\n", xl+1, yl, i);
	  }
	  
	} 
      }
    }
  }
  return n;
}

int mesh::up_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  int n = 0; 
  if (y != nyb*pow(2,level)) {
    int i;
    if (y % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->end(), found);
	  n++;
	  //found = H->search(xl-1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("UL.: (%d, %d): %d)\n", xl-1, yl, i);
	  //}
	  //found = H->search(xl+1, yl, i);
	  //if(found != NULL){
	  //  nb->insert(nb->begin(), found);
	    //printf ("UR.: (%d, %d): %d)\n", xl+1, yl, i);
	  //}
	  //cout << "FOUND in coarser" << endl;
	  //printf ("U.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (y % 2 == 0 || (y % 2 == 1 && i == -1)) {//não encontrou uma célula à esquerda mais grossa 
      for (i = level; i < number_of_levels; i++) {
	yl = (int) (pow (2, i - level) * y);
	yl = yl + pow(2, i - level); //procuramos pelos vizinhos à direita de c e a referência de c é o canto inferior esquerdo. Temos que "pular" as células nível atual (+ fino) dentro do próprio c 
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  xl = (int) (pow (2, i - level) * x) + j;
	  //printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->end(), found);
	    n++;
	    //found = H->search (xl-1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf("UL.: (%d, %d): %d)\n", xl-1, yl, i);
	    //}
	    //found = H->search (xl+1, yl, i);
	    //if (found != NULL) {
	    //  nb->insert (nb->begin(), found);
	      //printf("UR.: (%d, %d): %d)\n", xl+1, yl, i);
	    //}
	    //cout << "FOUND in finer" << endl;
	    //printf ("U.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
  return n;
}

list <cell *> * mesh::neighbours (cell * c) {
  list <cell *> * nb = new list <cell*>;
  int nw = this->left_neighbours(c, nb);
  c->set_cell_nvw(nw);
  int ns = this->down_neighbours(c, nb);
  c->set_cell_nvs(ns);
  int ne = this->right_neighbours(c, nb);
  c->set_cell_nve(ne);
  int nn = this->up_neighbours(c, nb);
  c->set_cell_nvn(nn);
  nb->insert(nb->begin(), c);
  return nb;
}

list <cell *> * mesh::neighbours_e(cell * c) {
  list <cell *> * nb = new list <cell*>;
  this->right_neighbours(c, nb);
  return nb;
}

list <cell *> * mesh::neighbours_w(cell * c) {
  list <cell *> * nb = new list <cell*>;
  this->left_neighbours(c, nb);
  return nb;
}

list <cell *> * mesh::neighbours_s(cell * c) {
  list <cell *> * nb = new list <cell*>;
  this->down_neighbours(c, nb);
  return nb;
}

list <cell *> * mesh::neighbours_n(cell * c) {
  list <cell *> * nb = new list <cell*>;
  this->up_neighbours(c, nb);
  return nb;
}


/* insere vizinhas na ordem */

void mesh::insert_neighbour(cell * c, list <cell *> * nb){
  unsigned int ct = 0;
  for(list <cell *>::iterator it = nb->begin(); it !=nb->end(); it++){ 
    if(c->get_cell_index() < (*it)->get_cell_index()){
      nb->insert(it, c);
      break;
    }
    ct++;
  }
  if(ct == nb->size())
    nb->insert(nb->end(), c);
}


/* Células da fronteira */
list <cell *> * mesh::boundary_cells_south(int tbc) {
  int i, j;
  list <cell *> * bcs = new list <cell*>;
  for (int k = 0; k < number_of_levels; k++) {
    for (list <cell *>::iterator it = l->at(k)->begin(); it != l->at(k)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      cell * c = search(i, j, k);
     
      if (c != NULL){
	if(j == 0){
	  //c->print_cell();
	  bcs->insert(bcs->begin(), c);
	  c->set_cell_bc(tbc);
	}
      }
    }
  }
  return bcs;
}

list <cell *> * mesh::boundary_cells_north(int tbc) {
  int i, j;
  
  list <cell *> * bcn = new list <cell*>;
  for (int k = 0; k < number_of_levels; k++) {
    for (list <cell *>::iterator it = l->at(k)->begin(); it != l->at(k)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      cell * c = search(i, j, k);
      if (c != NULL){
	if(j == pow(2,k)*nyb - 1){
	  //c->print_cell();
	  bcn->insert(bcn->begin(), c);
	  c->set_cell_bc(tbc);
	}
      }
    }
  }
  return bcn;
}

list <cell *> * mesh::boundary_cells_east(int tbc) {
  int i, j;
  list <cell *> * bce = new list <cell*>;
  for (int k = 0; k < number_of_levels; k++) {
    for (list <cell *>::iterator it = l->at(k)->begin(); it != l->at(k)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      cell * c = search(i, j, k);
      if (c != NULL){
	if(i == 0){
	  bce->insert(bce->begin(), c);
	  c->set_cell_bc(tbc);
	}
      }
    }
  }
  return bce;
}

list <cell *> * mesh::boundary_cells_west(int tbc) {
  int i, j;
  
  list <cell *> * bcw = new list <cell*>;
  for (int k = 0; k < number_of_levels; k++) {
    for (list <cell *>::iterator it = l->at(k)->begin(); it != l->at(k)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      cell * c = search(i, j, k);
      if (c != NULL){
	if(i == pow(2,k)*nxb - 1){
	  bcw->insert(bcw->begin(), c);
	  c->set_cell_bc(tbc);
	}
      }
    }
  }
  return bcw;
}

void mesh::print_cell_list(list <cell *> * bc){
  for(list <cell *>::iterator it = bc->begin(); it !=bc->end(); it++) 
    (*it)->print_cell();
   
}

void mesh::print_weight_list(list <weight> bc){
  for(list <weight>::iterator it = bc.begin(); it !=bc.end(); it++) 
    printf("%d %f\n", (it)->get_weight_index(), (it)->get_weight_w());
   
}


/* Finite Difference Neighbours */

void mesh::left_neighbours_fd (cell * c, list <cell *> * nb) {
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
    
  if (i != 0) {
    /* Same level ! */
    cell * found = H->search(i-1,j,level);
    if(found != NULL){
      insert_neighbour(found, nb);
    }
    else{
    /* Finer level! */
      found = H->search (2*(i-1), 2*j, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i-1), 2*j+1, level+1);
      if (found != NULL){
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i-1)+1, 2*j, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i-1)+1, 2*j+1, level+1);
      if (found != NULL) {
	insert_neighbour(found, nb);
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2-1, level-1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search ((i+1)/2-1, (j+1)/2, level-1);
      if (found != NULL){
	insert_neighbour(found, nb);
      }
    }
  }
}

void mesh::right_neighbours_fd (cell * c, list <cell *> * nb) {
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  
  if (i != pow(2,level)*nxb - 1) {
    cell * found = H->search (i+1, j, level);
    if (found != NULL) {
      insert_neighbour(found, nb);
    }
    else{
      /* Finer level! */
      found = H->search (2*(i+1), 2*j, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i+1), 2*j+1, level+1);
      if (found != NULL){
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i+1)+1, 2*j+1, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search (2*(i+1)+1, 2*j, level+1);
      if (found != NULL) {
	insert_neighbour(found, nb);
      }
      /* Coarser */
      found = H->search ((i+1)/2, (j+1)/2-1, level-1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
      }
      found = H->search ((i+1)/2, (j+1)/2, level-1);
      if (found != NULL){
	insert_neighbour(found, nb);
      }
    }
  }
}

void mesh::down_neighbours_fd (cell * c, list <cell *> * nb) {
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  list <cell *> * nbd = new list <cell*>;
  
  if (j != 0) {
     /* Same level ! */
    cell * found = H->search(i,j-1,level);
    if(found != NULL){
      insert_neighbour(found, nbd);
      insert_neighbour(found, nb);
    }
    else{
      /* Finer level! */
      found = H->search (2*i, 2*(j-1), level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
      found = H->search (2*i+1, 2*(j-1), level+1);
      if (found != NULL){
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
      found = H->search (2*i, 2*(j-1)+1, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
      found = H->search (2*i+1, 2*(j-1)+1, level+1);
      if (found != NULL) {
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2-1, level-1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
      found = H->search ((i+1)/2, (j+1)/2-1, level-1);
      if (found != NULL){
	insert_neighbour(found, nb);
	insert_neighbour(found, nbd);
      }
    }
  }
  
  for(list <cell *>::iterator it = nbd->begin(); it != nbd->end(); it++){
    this->left_neighbours_fd(*it, nb);
    this->right_neighbours_fd(*it, nb);
  }

}

void mesh::up_neighbours_fd (cell * c, list <cell *> * nb) {
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  list <cell *> * nbu = new list <cell*>;
  if (j != pow(2,level)*nyb - 1) {
    cell * found = H->search (i, j+1, level);
    if (found != NULL) {
      insert_neighbour(found, nb);
      insert_neighbour(found, nbu);
    }
    else{
      /* Finer level! */
      found = H->search (2*i, 2*(j+1), level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
      found = H->search (2*i+1, 2*(j+1), level+1);
      if (found != NULL){
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
      found = H->search (2*i+1, 2*(j+1)+1, level+1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
      found = H->search (2*i, 2*(j+1)+1, level+1);
      if (found != NULL) {
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2, level-1);
      if (found != NULL){ 
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
      found = H->search ((i+1)/2, (j+1)/2, level-1);
      if (found != NULL){
	insert_neighbour(found, nb);
	insert_neighbour(found, nbu);
      }
    }
  }
  for(list <cell *>::iterator it = nbu->begin(); it != nbu->end(); it++){
    this->left_neighbours_fd(*it, nb);
    this->right_neighbours_fd(*it, nb);
  }
}

list <cell *> * mesh::neighbours_fd (cell * c) {
  list <cell *> * nb = new list <cell*>;
  nb->insert(nb->begin(), c);  
  this->left_neighbours_fd(c, nb);
  this->right_neighbours_fd(c, nb);
  this->down_neighbours_fd(c, nb);
  this->up_neighbours_fd(c, nb);
  
  return nb;
}

/* Siblings */

vector <cell *> * mesh:: siblings(cell * c){
  vector <cell *> * L = new vector <cell *>;

  //cie => canto inferior esquerdo
  cell * cie;
  if (c->get_cell_x() % 2 == 1 && c->get_cell_y() % 2 == 1){
    //c is a cell top right
    cie = search(c->get_cell_x() - 1, c->get_cell_y() - 1, c->get_cell_level());
  }
  else if (c->get_cell_x() % 2 == 0 && c->get_cell_y() % 2 == 1){
    //c is a cell top left
    cie = search(c->get_cell_x(), c->get_cell_y() - 1, c->get_cell_level());
  }
  else if (c->get_cell_x() % 2 == 1 && c->get_cell_y() % 2 == 0){
    //c is a cell bottom right
    cie = search(c->get_cell_x() - 1, c->get_cell_y(), c->get_cell_level());
  }
  else {
    cie = c;
  }

  if (cie != NULL){// && cie->get_cell_refine() != number_of_levels) {
    cell * c_sibling;
    
    //célula irmã do canto inferior direito
    c_sibling = search(cie->get_cell_x() + 1, cie->get_cell_y(), cie->get_cell_level());
    if (c_sibling != NULL)// && c_sibling->get_cell_refine() != number_of_levels)
      L->push_back(c_sibling);

    //célula irmã do canto superior esquerdo
    c_sibling = search(cie->get_cell_x(), cie->get_cell_y() + 1, cie->get_cell_level());
    if (c_sibling != NULL)// && c_sibling->get_cell_refine() != number_of_levels)
      L->push_back(c_sibling);
    
    //célula irmã do canto superior direito
    c_sibling = search(cie->get_cell_x() + 1, cie->get_cell_y() + 1, cie->get_cell_level());
    if (c_sibling != NULL)// && c_sibling->get_cell_refine() != number_of_levels)
      L->push_back(c_sibling);
  }
  
  //feito no fim, depois da inclusão das células irmãs
  L->push_back(cie);

  //foram encontrados 4 células irmãs
  if (L->size() == 4)
    return L;
  else {
    L->clear();
    return NULL;
  }
}

bool equal(cell * c1, cell * c2){
  if (c1 == c2)
    return true;
  return false;
}

list <cell*>::iterator mesh::merge(vector <cell *> * V, cell *c){
  list <cell *>::iterator ret;
    
  cell * new_c = new cell ((V->back())->get_cell_x() / 2, (V->back())->get_cell_y() / 2, (V->back())->get_cell_level() - 1);
  //new_c->set_cell_with_particle(0);
  insert(new_c);
  for (vector <cell *>::iterator it = V->begin(); it != V->end(); it++) {
    if (!equal(*it, c)){
      remove(*it);
    }
  }
  ret = remove(c);
  return ret;
}

/* Conta numero de celulas da malha */

int mesh::counting_mesh_cell(){
  int ncell;
  ncell = 0;
    
  for (int i = 0; i < number_of_levels; i++) {
    for (list <cell *>::reverse_iterator it = l->at(i)->rbegin(); it != l->at(i)->rend(); it++) {
      (*it)->set_cell_index(ncell);      
      //(*it)->print_cell();
      ncell++;
    }
  }
  return(ncell);
}

/* Identify neighbours */

int mesh::search_left_neighbours(cell * c) {
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  int lv = -1;
    
  if (i != 0) {
    /* Same level ! */
    cell * found = H->search(i-1,j,level);
    if(found != NULL){
      lv = level;
    }
    else{
    /* Finer level! */
      found = H->search (2*(i-1), 2*j, level+1);
      if (found != NULL){ 
	lv = level + 1;
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2-1, level-1);
      if (found != NULL){ 
	lv = level - 1;
      }
    }
  }
  return(lv);
}

int mesh::search_down_neighbours(cell * c){
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  int lv = -1;
  
  if (j != 0) {
     /* Same level ! */
    cell * found = H->search(i,j-1,level);
    if(found != NULL){
      lv = level;
    }
    else{
      /* Finer level! */
      found = H->search (2*i, 2*(j-1), level+1);
      if (found != NULL){ 
	lv = level + 1;
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2-1, level-1);
      if (found != NULL){ 
	lv = level - 1;
      }
    }
  }
  return(lv);
}

int mesh::search_right_neighbours(cell * c){
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  int lv = -1;
  
  if (i != pow(2,level)*nxb - 1) {
    cell * found = H->search (i+1, j, level);
    if (found != NULL) {
      lv = level;
    }
    else{
      /* Finer level! */
      found = H->search (2*(i+1), 2*j, level+1);
      if (found != NULL){ 
	lv = level + 1;
      }
      /* Coarser */
      found = H->search ((i+1)/2, (j+1)/2-1, level-1);
      if (found != NULL){ 
	lv = level - 1;
      }
    }
  }
  return(lv);
}

int mesh::search_up_neighbours(cell * c){
  int i = c->get_cell_x();
  int j = c->get_cell_y();
  int level = c->get_cell_level();
  int lv = -1;
  
  if (j != pow(2,level)*nyb - 1) {
    cell * found = H->search (i, j+1, level);
    if (found != NULL) {
      lv = level;
    }
    else{
      /* Finer level! */
      found = H->search (2*i, 2*(j+1), level+1);
      if (found != NULL){ 
	lv = level + 1;
      }
      /* Coarser level */
      found = H->search ((i+1)/2-1, (j+1)/2, level-1);
      if (found != NULL){ 
	lv = level - 1;
      }
    }
  }
  return(lv);
}

/* Triangular mesh using hash table */

void mesh::triangular_mesh(int * ncellp, int * nelemp, int * narep, int * nvertp){
  int ncell, i, j, id, ind1, ind2, level, vert1, vert2, vert3, vert4, afacex0, afacexn, afacey0, afaceyn, nelem, nare, nvert, tam;
  nelem = * nelemp;
  nvert = * nvertp;
  nare = * narep;
  ncell = * ncellp;
  cell * c;
  double xbegin, xend, ybegin, yend, dx, dy;
  vector<cell *> * vert;
  vector<float> *vx, *vy;
  vector<int> * av1, * av2, * bc;
  vector<int> * ev1, * ev2, *ev3;
  vector<cell *> * elem;
  FILE * wFile;
  wFile = fopen("tmesh.dat","w");
  /* Boundary condition: Dirichlet 1, Neumann 2, Fluxo+ 4, Fluxo- 3, Interior 0 */
  afacex0 = 3;
  afacexn = 4;
  afacey0 = 2;
  afaceyn = 2;
  
  //nvert = 0;
  tam = 0;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();

  vert = new vector<cell *>;
  av1 = new vector<int>;
  av2 = new vector<int>;
  bc = new vector<int>;
  vx = new vector<float>;
  vy = new vector<float>;
  
  //fprintf(wFile, " Vertices \n");
  for (int k = 0; k < number_of_levels; k++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, k));
    dy = fabs(yend - ybegin) / (nyb * pow(2, k));
    for (list <cell *>::reverse_iterator it = l->at(k)->rbegin(); it != l->at(k)->rend(); it++) {
      (*it)->set_cell_index(ncell);
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      vert->push_back(*it);
      vx->push_back(xbegin+i*dx);
      vy->push_back(ybegin+j*dy);
      if(i == pow(2,level)*nxb - 1){
	nvert++;
	/* vertice a esquerda */
	vert->push_back(*it);
	vx->push_back(xbegin+(i+1)*dx);
	vy->push_back(ybegin+j*dy);
      }
      if(j == pow(2,level)*nyb - 1){
	nvert++;
	/* vertice acima */
	vert->push_back(*it);
	vx->push_back(xbegin+i*dx);
	vy->push_back(ybegin+(j+1)*dy);
      }
      if((j == pow(2,level)*nyb - 1) && (i == pow(2,level)*nxb - 1)){
	nvert++;
	/* vertice canto superior esquerdo */
	vert->push_back(*it);
	vx->push_back(xbegin+(i+1)*dx);
	vy->push_back(ybegin+(j+1)*dy);
      }
      nvert++;
      ncell++;
      tam++;
    }
  }
  
  nelem = 2*tam;
  cout << "# elem = " << nelem << endl;
  ev1 = new vector<int> (nelem);
  ev2 = new vector<int> (nelem);
  ev3 = new vector<int> (nelem);
  elem = new vector<cell *> (tam);
  tam = 0;
  for (int k = 0; k < number_of_levels; k++) {
    for (list <cell *>::reverse_iterator it = l->at(k)->rbegin(); it != l->at(k)->rend(); it++) {
      
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      level = (*it)->get_cell_level();
      c = search(i, j, level);
      if (c != NULL){
	id = c->get_cell_index();
	elem->at(tam) = (*it);
	tam++;
      }
    }
  }
  
  vert1 = vert2 = vert3 = vert4 = -1;
  //nare = 0;
  for (int k = 0; k < tam; k++) {
    ind1 = 2*k;
    ind2 = ind1 + 1;
    c = elem->at(k);
    i = c->get_cell_x();
    j = c->get_cell_y();
    id = c->get_cell_index();
    level = c->get_cell_level();
    if(c != NULL){
      cout << "Elemento = " << k << endl;
      vert1 = search_vertice(i, j, nvert, id, level, vert);
      if(j != nyb - 1){
	vert2 = search_vertice(i+1, j, nvert, id, level, vert);
	vert3 = search_vertice(i, j+1, nvert, id, level, vert);
	vert4 = search_vertice(i+1, j+1, nvert, id, level, vert);
	/* Primeiro elemento */
	ev1->at(ind1) = vert1;
	ev2->at(ind1) = vert2;
	ev3->at(ind1) = vert3;
	cout << ind1 << " " << vert1 << " " << vert2 << " " << vert3 << endl;
	cout << ind2 << " " << vert2 << " " << vert3 << " " << vert4 << endl;
	/* segundo elemento */
	ev1->at(ind2) = vert2;
	ev2->at(ind2) = vert3;
	ev3->at(ind2) = vert4;
	
	/* primeira aresta */
	av1->push_back(vert1);
	av2->push_back(vert2);
	if(j == 0)
	  bc->push_back(afacey0);
	else
	  bc->push_back(0);
	cout << nare << " " << vert1 << " " << vert2 << endl;
	nare++;
	/* segunda aresta */
	av1->push_back(vert1);
	av2->push_back(vert3);
	if(i == 0)
	  bc->push_back(afacex0);
	else
	  bc->push_back(0);
	cout << nare << " " << vert1 << " " << vert3 << endl;
	nare++;
	/* terceira aresta */
	av1->push_back(vert2);
	av2->push_back(vert3);
	bc->push_back(0);
	cout << nare << " " << vert2 << " " << vert3 << endl;
	nare++;
	if(i == nxb - 1){
	  /* quarta aresta */
	  av1->push_back(vert2);
	  av2->push_back(vert4);
	  bc->push_back(afacexn);
	  cout << nare << " " << vert2 << " " << vert4 << endl;
	  nare++;
	}
      }
      else{
	if(i != nxb - 1){
	  vert2 = search_vertice(i+1, j, nvert, id, level, vert); // vert 3
	  vert3 = vert1 + 1; //search_vertice(i, j + 1, nvert, id, level, vert);
	  vert4 = vert2 + 1; //search_vertice(i + 1, j + 1, nvert, id, level, 
	  
	  // Primeiro elemento 
	  ev1->at(ind1) = vert1;
	  ev2->at(ind1) = vert2;
	  ev3->at(ind1) = vert3;
	  cout << ind1 << " " << vert1 << " " << vert2 << " " << vert3 << endl;
	  cout << ind2 << " " << vert2 << " " << vert3 << " " << vert4 << endl;
	  // segundo elemento 
	  ev1->at(ind2) = vert2;
	  ev2->at(ind2) = vert3;
	  ev3->at(ind2) = vert4;
	  
	  // primeira aresta 
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  bc->push_back(0);
	  cout << nare << " " << vert1 << " " << vert2 << endl;
	  nare++;
	  // segunda aresta 
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  cout << nare << " " << vert1 << " " << vert3 << endl;
	  nare++;
	  // terceira aresta 
	  av1->push_back(vert2);
	  av2->push_back(vert3);
	  bc->push_back(0);
	  cout << nare << " " << vert2 << " " << vert3 << endl;
	  nare++;
	  // quarta aresta 
	  av1->push_back(vert3);
	  av2->push_back(vert4);
	  bc->push_back(afaceyn);
	  cout << nare << " " << vert4 << " " << vert3 << endl;
	  nare++;
	}
	else{
	  vert2 = vert1 + 2; //search_vertice(i+1, j, nvert, id, level, vert); 
	  vert3 = vert1 + 1; //search_vertice(i, j + 1, nvert, id, level, vert)
	  vert4 = vert2 + 1; //search_vertice(i + 1, j + 1, nvert, id, level,
	  
	  // Primeiro elemento 
	  ev1->at(ind1) = vert1;
	  ev2->at(ind1) = vert2;
	  ev3->at(ind1) = vert3;
	  cout << ind1 << " " << vert1 << " " << vert2 << " " << vert3 << endl;
	  cout << ind2 << " " << vert2 << " " << vert3 << " " << vert4 << endl;
	  // segundo elemento 
	  ev1->at(ind2) = vert2;
	  ev2->at(ind2) = vert3;
	  ev3->at(ind2) = vert4;
	  
	  // primeira aresta 
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  bc->push_back(0);
	  cout << nare << " " << vert1 << " " << vert2 << endl;
	  nare++;
	  // segunda aresta 
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  bc->push_back(0);
	  cout << nare << " " << vert1 << " " << vert3 << endl;
	  nare++;
	  // terceira aresta 
	  av1->push_back(vert2);
	  av2->push_back(vert3);
	  bc->push_back(0);
	  cout << nare << " " << vert2 << " " << vert3 << endl;
	  nare++;
	  //quarta aresta
	  av1->push_back(vert3);
	  av2->push_back(vert4);
	  bc->push_back(3);
	  nare++;
	  // quinta aresta
	  av1->push_back(vert4);
	  av2->push_back(vert2);
	  bc->push_back(2);
	  nare++;
	}
      }
    }
  }
  fprintf(wFile,"%d %d %d\n", nelem, nvert, nare);
  cout << ncell << " " << nare << " " << nvert << " " << nelem << endl;
  *nelemp = nelem;
  *narep = nare;
  *nvertp = nvert;
  *ncellp = ncell;

  for(i = 0; i < nvert; i++)
    fprintf(wFile,"%d %f %f\n", i, vx->at(i), vy->at(i));
      
  for(i = 0; i < nelem; i++)
    fprintf(wFile, "%d %d %d\n", ev1->at(i), ev2->at(i), ev3->at(i));
  
  for(i = 0; i < nare; i++)
    fprintf(wFile, "%d %d %d\n", av1->at(i), av2->at(i), bc->at(i));

  /* Aqui imprimir elem regiao -1 */
  for(i = 0; i < nelem; i++)
    fprintf(wFile, "%d %d %d\n", i, -1, -1);
    
  delete vert;
  delete ev1;
  delete ev2;
  delete ev3;
  delete av1;
  delete av2;
  delete elem;
  fclose(wFile);  
}

void mesh::triangular_mesh_refined(int * ncellp, int * nelemp, int * narep, int * nvertp){
  int ncell, i, j, vert1, vert2, vert3, vert4, vert5, vert6, vert7, lvu, lvd;
  int afacex0, afacexn, afacey0, afaceyn, nelem, nare, nvert, tam, lvl, lvr;
  float xv, yv;
  nelem = * nelemp;
  nvert = * nvertp;
  nare = * narep;
  ncell = * ncellp;
  cell * c;
  double xbegin, xend, ybegin, yend, dx, dy;
  vector<cell *> * vert;
  vector<float> *vx, *vy;
  vector<int> * av1, * av2, * bc;
  vector<int> * ev1, * ev2, *ev3, *indv;
  vector<cell *> * elem;
  FILE * wFile;
  wFile = fopen("tmeshr.dat","w");
  /* Boundary condition: Dirichlet 1, Neumann 2, Fluxo+ 4, Fluxo- 3, Interior 0 */
  afacex0 = 3;
  afacexn = 4;
  afacey0 = 2;
  afaceyn = 2;
  
  //nvert = 0;
  tam = 0;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();

  indv = new vector<int>;
  vert = new vector<cell *>;
  av1 = new vector<int>;
  av2 = new vector<int>;
  ev1 = new vector<int>;
  ev2 = new vector<int>;
  ev3 = new vector<int>;
  bc = new vector<int>;
  vx = new vector<float>;
  vy = new vector<float>;
  elem = new vector<cell *>;
  
  for (int k = 0; k < number_of_levels; k++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, k));
    dy = fabs(yend - ybegin) / (nyb * pow(2, k));
    for (list <cell *>::reverse_iterator it = l->at(k)->rbegin(); it != l->at(k)->rend(); it++) {
      (*it)->set_cell_index(ncell);
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      vx->push_back(xbegin+i*dx);
      vy->push_back(ybegin+j*dy);
      nelem += 2;
      indv->push_back(nvert);
      c = search(i, j, k);
      if (c != NULL){
	elem->push_back(*it);
	lvl = search_left_neighbours(c);
	lvr = search_right_neighbours(c);
	lvd = search_down_neighbours(c);
	lvu = search_up_neighbours(c);
	if(lvl > k){
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+i*dx);
	  vy->push_back(ybegin+(j+0.5)*dy);
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+(i+0.5)*dx);
	  vy->push_back(ybegin+(j+0.5)*dy);
	  nelem += 3;
	  if(lvd > k){
	    nvert++;
	    indv->push_back(nvert);
	    vx->push_back(xbegin+(i+0.5)*dx);
	    vy->push_back(ybegin+j*dy);
	    nelem++;
	  }
	  if(lvr > k)
	    nelem++;
	  if(lvu > k)
	    nelem++;
	}
	else if(lvr > k){
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+(i+0.5)*dx);
	  vy->push_back(ybegin+(j+0.5)*dy);
	  nelem += 3;
	  if(lvu > k)
	    nelem++;
	  if(lvd > k){
	    nvert++;
	    indv->push_back(nvert);
	    vx->push_back(xbegin+(i+0.5)*dx);
	    vy->push_back(ybegin+j*dy);
	    nelem++;
	  }
	}
	else if(lvd > k){
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+(i+0.5)*dx);
	  vy->push_back(ybegin+j*dy);
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+(i+0.5)*dx);
	  vy->push_back(ybegin+(j+0.5)*dy);
	  nelem += 3;
	  if(lvu > k)
	    nelem++;
	}
	else if(lvu > k){
	  nvert++;
	  indv->push_back(nvert);
	  vx->push_back(xbegin+(i+0.5)*dx);
	  vy->push_back(ybegin+(j+0.5)*dy);
	  nelem += 3;
	}
      }
      if(i == pow(2,k)*nxb - 1){
	nvert++;
	indv->push_back(nvert);
	//vertice a esquerda 
	vx->push_back(xbegin+(i+1)*dx);
	vy->push_back(ybegin+j*dy);
      }
      if(j == pow(2,k)*nyb - 1){
	nvert++;
	indv->push_back(nvert);
	// vertice acima 
	vx->push_back(xbegin+i*dx);
	vy->push_back(ybegin+(j+1)*dy);
      }
      if((j == pow(2,k)*nyb - 1) && (i == pow(2,k)*nxb - 1)){
	nvert++;
	indv->push_back(nvert);
	// vertice canto superior esquerdo 
       	vx->push_back(xbegin+(i+1)*dx);
	vy->push_back(ybegin+(j+1)*dy);
      }
      nvert++;
      ncell++;
      tam++;
    }
  } 
  cout << "# elem = " << nelem << " #nvert = " << nvert << " #cell = " << tam << endl;
  ordena_vertice(nvert, indv, vx, vy);
  
  vert1 = vert2 = vert3 = vert4 = -1;
  vert5 = vert6 = vert7 = -1;
  
  for (int k = 0; k < number_of_levels; k++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, k));
    dy = fabs(yend - ybegin) / (nyb * pow(2, k));
    for (list <cell *>::reverse_iterator it = l->at(k)->rbegin(); it != l->at(k)->rend(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      xv = xbegin + i*dx;
      yv = ybegin + j*dy;
      c = search(i, j, k);
      if (c != NULL){
        vert1 = search_vertice_index(xv, yv, nvert, vx, vy, indv);
	vert2 = search_vertice_index(xv + dx, yv, nvert, vx, vy, indv);
	vert3 = search_vertice_index(xv, yv + dy, nvert, vx, vy, indv);
	vert4 = search_vertice_index(xv + dx, yv + dy, nvert, vx, vy, indv);
	lvl = search_left_neighbours(c);
	lvr = search_right_neighbours(c);
	lvd = search_down_neighbours(c);
	lvu = search_up_neighbours(c);
	if(lvl > k){
	  vert6 = search_vertice_index(xv, yv + 0.5*dy, nvert, vx, vy, indv);
	  vert5 = search_vertice_index(xv + 0.5*dx, yv + 0.5*dy, nvert, vx, vy, indv);
	  // elemento 1
	  ev1->push_back(vert1);
	  ev2->push_back(vert2);
	  ev3->push_back(vert5);
	  // elemento 2
	  ev1->push_back(vert2);
	  ev2->push_back(vert4);
	  ev3->push_back(vert5);
	  // elemento 3
	  ev1->push_back(vert4);
	  ev2->push_back(vert3);
	  ev3->push_back(vert5);
	  // elemento 4
	  ev1->push_back(vert3);
	  ev2->push_back(vert5);
	  ev3->push_back(vert6);
	  // elemento 5
	  ev1->push_back(vert1);
	  ev2->push_back(vert5);
	  ev3->push_back(vert6);
	  // aresta 1
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 2 
	  av1->push_back(vert1);
	  av2->push_back(vert6);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 3
	  av1->push_back(vert6);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 4
	  av1->push_back(vert1);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 5
	  av1->push_back(vert4);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 6
	  av1->push_back(vert2);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 7
	  av1->push_back(vert3);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 8
	  av1->push_back(vert6);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  if(i == pow(2,k)*nxb - 1){
	    // aresta 9
	    av1->push_back(vert2);
	    av2->push_back(vert4);
	    bc->push_back(afacexn);
	    nare++;
	  }
	  if(lvd > k){
	    vert7 = search_vertice_index(xv + 0.5*dx, yv + dy, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	    if(j == pow(2,k)*nyb - 1){
	      // aresta 4
	      av1->push_back(vert3);
	      av2->push_back(vert4);
	      bc->push_back(afaceyn);
	      nare++;
	    }
	  }
	  if(lvr > k){
	    vert7 = search_vertice_index(xv + dx, yv + 0.5*dy, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	  }
	  if(lvu > k){
	    vert7 = search_vertice_index(xv + 0.5*dx, yv, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	  }
	}
	else if(lvr > k){
	  vert6 = search_vertice_index(xv + dx, yv + 0.5*dy, nvert, vx, vy, indv);
	  vert5 = search_vertice_index(xv + 0.5*dx, yv + 0.5*dy, nvert, vx, vy, indv);
	  // elemento 1
	  ev1->push_back(vert1);
	  ev2->push_back(vert2);
	  ev3->push_back(vert5);
	  // elemento 2
	  ev1->push_back(vert2);
	  ev2->push_back(vert6);
	  ev3->push_back(vert5);
	  // elemento 3
	  ev1->push_back(vert4);
	  ev2->push_back(vert6);
	  ev3->push_back(vert5);
	  // elemento 4
	  ev1->push_back(vert3);
	  ev2->push_back(vert5);
	  ev3->push_back(vert4);
	  // elemento 5
	  ev1->push_back(vert1);
	  ev2->push_back(vert5);
	  ev3->push_back(vert3);
	  // aresta 1
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 2 
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 3
	  av1->push_back(vert5);
	  av2->push_back(vert3);
	  bc->push_back(0);
	  nare++;
	  // aresta 4
	  av1->push_back(vert1);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 5
	  av1->push_back(vert4);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 6
	  av1->push_back(vert2);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 7
	  av1->push_back(vert6);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  if(lvd > k){
	    vert7 = search_vertice_index(xv + 0.5*dx, yv + dy, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	  }
	  if(lvu > k){
	    vert7 = search_vertice_index(xv + 0.5*dx, yv, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	  }
	}
	else if(lvd > k){
	  vert6 = search_vertice_index(xv + 0.5*dx, yv, nvert, vx, vy, indv);
	  vert5 = search_vertice_index(xv + 0.5*dx, yv + 0.5*dy, nvert, vx, vy, indv);
	  // elemento 1
	  ev1->push_back(vert1);
	  ev2->push_back(vert6);
	  ev3->push_back(vert5);
	  // elemento 2
	  ev1->push_back(vert2);
	  ev2->push_back(vert6);
	  ev3->push_back(vert5);
	  // elemento 3
	  ev1->push_back(vert4);
	  ev2->push_back(vert2);
	  ev3->push_back(vert5);
	  // elemento 4
	  ev1->push_back(vert3);
	  ev2->push_back(vert5);
	  ev3->push_back(vert4);
	  // elemento 5
	  ev1->push_back(vert1);
	  ev2->push_back(vert5);
	  ev3->push_back(vert3);
	  // aresta 1
	  av1->push_back(vert1);
	  av2->push_back(vert6);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 2 
	  av1->push_back(vert2);
	  av2->push_back(vert6);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 3
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 4
	  av1->push_back(vert1);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 5
	  av1->push_back(vert4);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 6
	  av1->push_back(vert2);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 7
	  av1->push_back(vert3);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 8
	  av1->push_back(vert6);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  if(j == pow(2,k)*nyb - 1){
	    // aresta 9
	    av1->push_back(vert3);
	    av2->push_back(vert4);
	    bc->push_back(afaceyn);
	    nare++;
	  }
	  if(lvu > k){
	    vert7 = search_vertice_index(xv + 0.5*dx, yv, nvert, vx, vy, indv);
	    av1->push_back(vert5);
	    av2->push_back(vert7);
	    bc->push_back(0);
	    nare++;
	  }
	}
	else if(lvu > k){
	  vert6 = search_vertice_index(xv + 0.5*dx, yv + dy, nvert, vx, vy, indv);
	  vert5 = search_vertice_index(xv + 0.5*dx, yv + 0.5*dy, nvert, vx, vy, indv);
	  // elemento 1
	  ev1->push_back(vert1);
	  ev2->push_back(vert2);
	  ev3->push_back(vert5);
	  // elemento 2
	  ev1->push_back(vert1);
	  ev2->push_back(vert3);
	  ev3->push_back(vert5);
	  // elemento 3
	  ev1->push_back(vert4);
	  ev2->push_back(vert2);
	  ev3->push_back(vert5);
	  // elemento 4
	  ev1->push_back(vert3);
	  ev2->push_back(vert5);
	  ev3->push_back(vert6);
	  // elemento 5
	  ev1->push_back(vert4);
	  ev2->push_back(vert5);
	  ev3->push_back(vert6);
	  // aresta 1
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 2 
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 3
	  av1->push_back(vert5);
	  av2->push_back(vert3);
	  bc->push_back(0);
	  nare++;
	  // aresta 4
	  av1->push_back(vert1);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 5
	  av1->push_back(vert4);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 6
	  av1->push_back(vert2);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	  // aresta 7
	  av1->push_back(vert6);
	  av2->push_back(vert5);
	  bc->push_back(0);
	  nare++;
	}
	else{
	  // elemento 1
	  ev1->push_back(vert1);
	  ev2->push_back(vert2);
	  ev3->push_back(vert4);
	  // elemento 2
	  ev1->push_back(vert1);
	  ev2->push_back(vert3);
	  ev3->push_back(vert4);
	  // aresta 1
	  av1->push_back(vert1);
	  av2->push_back(vert2);
	  if(j == 0)
	    bc->push_back(afacey0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 2 
	  av1->push_back(vert1);
	  av2->push_back(vert3);
	  if(i == 0)
	    bc->push_back(afacex0);
	  else
	    bc->push_back(0);
	  nare++;
	  // aresta 3 
	  av1->push_back(vert1);
	  av2->push_back(vert4);
	  bc->push_back(0);
	  nare++;
	  if(i == pow(2,k)*nxb - 1){
	    // aresta 4
	    av1->push_back(vert2);
	    av2->push_back(vert4);
	    bc->push_back(afacexn);
	    nare++;
	  }
	  if(j == pow(2,k)*nyb - 1){
	    // aresta 4
	    av1->push_back(vert3);
	    av2->push_back(vert4);
	    bc->push_back(afaceyn);
	    nare++;
	  }
	}
      }
    }
  }
  fprintf(wFile,"%d %d %d\n", nelem, nvert, nare);
  cout << ncell << " " << nelem << " " << nvert << " " << nare << endl;
 
  for(i = 0; i < nvert; i++)
    fprintf(wFile,"%d %f %f\n", indv->at(i), vx->at(i), vy->at(i));
 
  
  for(i = 0; i < nelem; i++)
    fprintf(wFile, "%d %d %d\n", ev1->at(i), ev2->at(i), ev3->at(i));

  for(i = 0; i < nare; i++)
    fprintf(wFile, "%d %d %d\n", av1->at(i), av2->at(i), bc->at(i));

  // Aqui imprimir elem regiao -1 
  //for(i = 0; i < nelem; i++){
  //  fprintf(wFile, "%d %d %d\n", i, 0, -1);
  //}
  delete vert;
  delete ev1;
  delete ev2;
  delete ev3;
  delete av1;
  delete av2;
  delete elem;
  fclose(wFile);  
}

/* Busca vertice na lista de vertices */

int mesh::search_vertice(int i, int j, int n, int id, int level, vector <cell *> * vert){
  int ivert = -1, k, ii, jj, l, m;
  
  k = 0;
  while(k < n && vert->at(k)->get_cell_index() < id)
    k++;

  ii = vert->at(k)->get_cell_x();
  jj = vert->at(k)->get_cell_y();
    
  if((i == ii) && (j == jj))
    ivert = k;
  if((i == ii) && (j == jj + 1)){
    l = k;
    while(l < n && vert->at(l)->get_cell_y() < j)
      l++;
    ivert = l + i; // Tratar quando houver refinamento
    if((j == nyb - 1) && (i != 0))
      ivert = l + 2*i; 
  }
  if((i == ii + 1) && (j == jj)){
    l = k;
    while(l < n && vert->at(l)->get_cell_x() < i)
      l++;
    if(i != nxb)
      ivert = l;
    else{
      ivert = i + j*(nxb+1); // Tratar quando houver refinamento
    }
  }
  if((i == ii + 1) && (j == jj + 1)){
    l = k;
    while(l < n && vert->at(l)->get_cell_y() < j)
      l++;
    m = l + i;
    while(m < n && vert->at(m)->get_cell_x() < i)
      m++;
    
    if(i != nxb)
      ivert = m;
    else{
      ivert = i + j*(nxb+1);
      if(j == nyb - 1){
	ivert = i + j*(nxb+1) + nxb; //Tratar caso com refinamento
      }
    }

    
    
  }
   
  //cout << "vertice " << i << " " << j << " " << " " << ivert << endl;
  
  return(ivert);
}

/* Ordena indices dos vertices em relacao a vx e vy */

void mesh::ordena_vertice(int nvert, vector <int> * vind, vector <float> * vx, vector <float> * vy){
  int i, j, k, ind, ia = 0;
  float x, y;
  for(j = 1; j < nvert; j++){
    x = vx->at(j);
    y = vy->at(j);
    ind = vind->at(j);
    for(i = j - 1; i >= 0 && vx->at(i) > x; i--){
      vx->at(i+1) = vx->at(i);
      vy->at(i+1) = vy->at(i);
      vind->at(i+1) = vind->at(i);
    }
    vx->at(i+1) = x;
    vy->at(i+1) = y;
    vind->at(i+1) = ind;
  }
  
  for(j = 1; j < nvert; j++){
    if(vx->at(j-1) != vx->at(j)){
      for(k = ia + 1; k < j; k++){
	x = vx->at(k);
	y = vy->at(k);
	ind = vind->at(k);
	for(i = k - 1; i >= ia && vy->at(i) > y; i--){
	  vx->at(i+1) = vx->at(i);
	  vy->at(i+1) = vy->at(i);
	  vind->at(i+1) = vind->at(i);
	}
	vx->at(i+1) = x;
	vy->at(i+1) = y;
	vind->at(i+1) = ind;
      }
      ia = j;
     }
  }
    
}

/* Busca indice do vertice na lista de vertices ordenada */

int mesh::search_vertice_index(float x, float y, int nvert, vector <float> * vx, vector <float> * vy, vector <int> * vind){
  int ivert = -1, k, i;
  
  k = 0;
  while(k < nvert && vx->at(k) < x)
    k++;
  i = k;
  while(i < nvert && vy->at(i) < y)
    i++;
  ivert = vind->at(i);
  return(ivert);
} 

/* cria vetor com endereco da celulas na ordem da numeração */

void mesh::mesh_adress(int ncell, vector <cell*> * elem){
  
  int k, i, j, level, id = -1;
  cell * c = NULL;
  
  //cout << "ELEM " << elem->size() << endl;
  for (k = 0; k < number_of_levels; k++) {
    for (list <cell *>::iterator it = l->at(k)->begin(); it != l->at(k)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      level = (*it)->get_cell_level();
      c = search(i, j, level);
      if (c != NULL){
	id = c->get_cell_index();
	elem->at(id) = (*it);
	//cout << id << " " << elem[id]->get_cell_x() << "," << elem[id]->get_cell_y() << endl;
	
      }      
    }
  }
  
}

/*Vizinhas de todas as celulas da malha */

int mesh::neighbours_all_cell(){
  int i, j, level, idva, idv, idx;
  int ncellv = 0;
  cell * c;
  
  for (int ll = 0; ll < number_of_levels; ll++) {
    for (list <cell *>::iterator it = l->at(ll)->begin(); it != l->at(ll)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      level = (*it)->get_cell_level();
      idx = (*it)->get_cell_index();
      c = search(i, j, level);
      if (c != NULL) {
	list <cell *> * lnew = neighbours(c);
	printf("%d ",idx);
	//(*it)->print_cell();
	idva = -1;
	for (list <cell *>::iterator ite = lnew->begin(); ite != lnew->end(); ite++) {
	  idv = (*ite)->get_cell_index();
	  printf("%d ", idv);
	  //(*ite)->print_cell();
	  if(idv != idva)
	    ncellv++;
	  idva = idv;
	}
	lnew->clear();
	printf("\n");
      }
      else
	{
	  cout << "There is no such cell!" << endl;
	}
    }
  }
  return(ncellv);
}

/* Funcao definida na malha */

void mesh::calculation_function (double (* f) (double x, double y), vector <double> * fvalue){
  double xbegin, ybegin, xend, yend, dx, dy, xmiddle, ymiddle;
  int i, j;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();
  
  for (int ll = 0; ll < number_of_levels; ll++) {
    dx = fabs(xend - xbegin) / (nxb * pow(2, ll));
    dy = fabs(yend - ybegin) / (nyb * pow(2, ll));
    for (list <cell *>::iterator it = l->at(ll)->begin(); it != l->at(ll)->end(); it++) {
      i = (*it)->get_cell_x();
      j = (*it)->get_cell_y();
      
      xmiddle = xbegin + (i + 0.5) * dx;//(delta_0 / pow(2, i)))
      ymiddle = ybegin + (j + 0.5) * dy;//(delta_0 / pow(2, i)));
      //cout << xmiddle << " " << ymiddle << " " << endl;
      
      fvalue->push_back((*f) (xmiddle, ymiddle));
      
      //printf ("%d \t %.7f \t %.7f \t %.7f\n %", ncell, xmiddle, ymiddle, (*f)(xmiddle, ymiddle));
      
    }
    //printf ("\n");
  }
  
}

/* Inicializa RHS - Boundary conditions */

void mesh::rhs_dirichlet_boundary_conditions (double (*f) (double x, double y), vector<double> * fvalue) {
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

      //      cout << ll << " " << i << " " << j << " " << xmiddle << " " << ymiddle << " " << xbegin << " " << xend << " " << dx << endl;
      
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
} 
 

/* Norma do maximo de um vetor */

double mesh::normmax(vector <double> * x, int ncell){
  int i, aux, rmax = 0.0;
  for(i = 0; i < ncell; i++){
    aux = fabs(x->at(i));
    if(aux > rmax)
      rmax = aux;
  }
  return(rmax);
}

/* Norma 2 do erro */

double mesh::erron2(vector <double> * x, vector <double> * xa, int ncell){
  int i, soma = 0.0, erro_norm2;
  for(i = 0; i < ncell; i++)
    soma += (x->at(i) - xa->at(i))*(x->at(i) - xa->at(i));
  erro_norm2 = sqrt(soma);
  return(erro_norm2);
}

/* imprime coeficientes nao nulos da matriz A */

void mesh::print_matrix(vector <int> * IA, vector <int> * JA, vector <double> * A, int ncell, int ncellv){
  int i;
  
  for (i = 0; i < ncell + 1; i++){
    //cout << "IA(" << i << ")= " << IA->at(i) << endl;
    cout << IA->at(i) << endl;
  }
  cout << endl;
  for (i = 0; i < ncellv; i++){
    //cout << "A(" << i << ")= " << A->at(i) << endl;
    cout << A->at(i) << endl;
  }
  cout << endl;
  for (i = 0; i < ncellv; i++){
    //cout << "JA(" << i << ")= " << JA->at(i) << endl;
    cout << JA->at(i) << endl;
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vector <double> mesh::max_propriedades(){
  vector <double> pmax;
  double velu, velv, phi;
  pmax.push_back(0.0);
  pmax.push_back(0.0);
  pmax.push_back(0.0);
 
  for (int i = 0; i < number_of_levels; i++) {
    
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      velu = (*it)->get_cell_velu();
      velv = (*it)->get_cell_velv();
      phi = (*it)->get_cell_phi();
      if(velu > pmax[0])
	pmax[0] = velu;
      if(velv > pmax[1])
	pmax[1] = velv;
      if(phi > pmax[2])
	pmax[2] = phi;
    }
  }
  return pmax;
}
