#include <cmath>
#include <cassert>
#include "mesh.h"

mesh::mesh(){
  l = NULL;
  number_of_levels = -1;
  max_dimension_by_level = NULL;
}

mesh::mesh(dominio * D, int n_levels, vector <int> * max_dimension_by_level){
  this->D = D;
  //o tamanho da tabela hash vai depender do máximo número de células do nível mais fino, isto é:
  //(1/100)*max_dim_by_level[n_level - 1] => (x * 100) % de max_dim_by_level[n_level - 1]
  double x = 0.2;
  H = new hash_table((int)(x * max_dimension_by_level->at(n_levels - 1)));
  //H = new hash_table(1000000);
  l = new vector<list<cell *> * >;
  for (int i = 0; i < n_levels; i++){
    list<cell *> * l_by_level= new list <cell *>;
    l->push_back(l_by_level); 
  }
  number_of_levels = n_levels;
  this->max_dimension_by_level = max_dimension_by_level;
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

cell * mesh::search(int x, int y, int level) {
  return H->search(x, y, level);
}

void mesh::print_mesh() {
  for (int i = 0; i < number_of_levels; i++)
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++)
      cout << "(" << (*it)->get_cell_x() << ", " << (*it)->get_cell_y() << "): " << (*it)->get_cell_level() << endl;
  
}

hash_table * mesh::get_hash_table() {
  return H;
}


list <cell *>::iterator mesh::split_and_insert(cell * c) {
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
  /****************************************************************************/
  
  double_hash_table * internal_HT = new double_hash_table(4 * H->get_number_cell());
  
  int nzones = 0;
  for (int i = 0; i < number_of_levels; i++) {
    //printf ("Nível: %d\n", i);
    dx = fabs(xend - xbegin) / (32 * pow(2, i));
    dy = fabs(yend - ybegin) / (32 * pow(2, i));
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
      delete cmiddle;
      nzones++;
    }
  }
  cout << nzones << " " << ncells << endl;
  assert (nzones == ncells);

  /********************silo code*******************************/
  DBfile *dbfile = NULL;
  /* Open the Silo file */
  char palavra [100];
  for (int i = 0; i < 100; i++)
    palavra [0] = '\0';

  sprintf(palavra, "/home/alvaro/Desktop/data/basic%4.3lf.silo", 1 + tempo); 
  dbfile = DBCreate(palavra, DB_CLOBBER, DB_LOCAL,"Comment about the data", DB_HDF5);
  if(dbfile == NULL)
    {
      fprintf(stderr, "Could not create Silo file!\n");
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
  assert(DBPutUcdvar1 (dbfile, "cellmiddles", "meshBLA", cellmiddles, nzones, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL) == 0);
      
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

void mesh::calculation_exact_function (double (* f) (double_cell * c)){
  double xbegin, ybegin, xend, yend, dx, dy, xmiddle, ymiddle;
  double_cell * c;
  xbegin = get_dominio()->get_xbegin();
  ybegin = get_dominio()->get_ybegin();
  xend = get_dominio()->get_xend();
  yend = get_dominio()->get_yend();
  for (int i = 0; i < number_of_levels; i++) {
    dx = fabs(xend - xbegin) / (32 * pow(2, i));
    dy = fabs(yend - ybegin) / (32 * pow(2, i));
    for (list <cell *>::iterator it = l->at(i)->begin(); it != l->at(i)->end(); it++) {
      xmiddle = xbegin + (((*it)->get_cell_x() * dx) / 2.);//(delta_0 / pow(2, i)));
      ymiddle = ybegin + (((*it)->get_cell_y() * dy) / 2.);//(delta_0 / pow(2, i)));
      
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

void mesh::left_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  
  dominio * D = get_dominio();
  
  if (x != D->get_xbegin()) {
    int i;
    if (x % 2 == 0){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	xl = xl - 1;
	//printf ("(%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);
	  cout << "FOUND in coarser" << endl;
	  printf ("L.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (x % 2 == 1 || (x % 2 == 0 && i == -1)) {//não encontrou uma célula à esquerda mais grossa 
      for (i = level; i < number_of_levels; i++) {
	xl = (int) (pow (2, i - level) * x);
	xl = xl - 1;
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  yl = (int) (pow (2, i - level) * y) + j;
	  //printf ("(%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->begin(), found);
	    cout << "FOUND in finer" << endl;
	    printf ("L.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
}

void mesh::down_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  
  dominio * D = get_dominio();
  
  if (y != D->get_ybegin()) {
    int i;
    if (y % 2 == 0){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl - 1;
	//printf ("(%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);
	  cout << "FOUND in coarser" << endl;
	  printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	  break;
	}
      }
    }
    if (y % 2 == 1 || (y % 2 == 0 && i == -1)) {//não encontrou uma célula à esquerda mais grossa 
      for (i = level; i < number_of_levels; i++) {
	yl = (int) (pow (2, i - level) * y);
	yl = yl - 1;
	for (int j = 0; j <= (int) pow (2, i - level) - 1; j++) {
	  xl = (int) (pow (2, i - level) * x) + j;
	  //printf ("(%d, %d): %d)\n", xl, yl, i);
	  cell * found = H->search (xl, yl, i);
	  if (found != NULL) {
	    nb->insert (nb->begin(), found);
	    cout << "FOUND in finer" << endl;
	    printf ("D.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
}

void mesh::right_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  
  dominio * D = get_dominio();
  
  if (x != D->get_xend()) {
    int i;
    if (x % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	xl = xl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);
	  cout << "FOUND in coarser" << endl;
	  printf ("R.: (%d, %d): %d)\n", xl, yl, i);
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
	    nb->insert (nb->begin(), found);
	    cout << "FOUND in finer" << endl;
	    printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
}

void mesh::up_neighbours (cell * c, list <cell *> * nb) {
  int x = c->get_cell_x();
  int y = c->get_cell_y();
  int level = c->get_cell_level();
  int xl, yl;
  
  dominio * D = get_dominio();
  
  if (y != D->get_yend()) {
    int i;
    if (y % 2 == 1){
      for (i = level - 1; i >= 0; i--){
	xl = (int) (x / pow (2, level - i));
	yl = (int) (y / pow (2, level - i));
	yl = yl + 1;
	//printf ("R.: (%d, %d): %d)\n", xl, yl, i);
	cell * found = H->search (xl, yl, i);
	if (found != NULL) {
	  nb->insert (nb->begin(), found);
	  cout << "FOUND in coarser" << endl;
	  printf ("U.: (%d, %d): %d)\n", xl, yl, i);
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
	    nb->insert (nb->begin(), found);
	    cout << "FOUND in finer" << endl;
	    printf ("U.: (%d, %d): %d)\n", xl, yl, i);
	  }
	} 
      }
    }
  }
}

list <cell *> * mesh::neighbours (cell * c) {
  list <cell *> * nb = new list <cell*>;

  cout << nb->size() << endl;
  cout << "left" << endl;
  this->left_neighbours(c, nb);
  cout << nb->size() << endl;
  cout << "right" << endl;
  this->right_neighbours(c, nb);
  cout << nb->size() << endl;
  cout << "down" << endl;
  this->down_neighbours(c, nb);
  cout << nb->size() << endl;
  cout << "up" << endl;
  this->up_neighbours(c, nb);

  return nb;
}

vector <cell *> * mesh::siblings (cell * c) {
  vector <cell *> * L = new vector <cell *>;

  //cie => canto inferior esquerdo
  cell * cie;
  if ((c->get_cell_kind())[0] == 't' && (c->get_cell_kind())[1] == 'r'){
    cie = search(c->get_cell_x() - 1, c->get_cell_y() - 1, c->get_cell_level());
  }
  else if ((c->get_cell_kind())[0] == 't' && (c->get_cell_kind())[1] == 'l'){
    cie = search(c->get_cell_x(), c->get_cell_y() - 1, c->get_cell_level());
  }
  else if ((c->get_cell_kind())[0] == 'b' && (c->get_cell_kind())[1] == 'r'){
    cie = search(c->get_cell_x() - 1, c->get_cell_y(), c->get_cell_level());
  }
  else {
    cie = c;
  }

  if (cie != NULL) {
    cell * c_sibling;
    
    //célula irmã do canto inferior direito
    c_sibling = search(cie->get_cell_x() + 1, cie->get_cell_y(), cie->get_cell_level());
    if (c_sibling != NULL)
      L->push_back(c_sibling);

    //célula irmã do canto superior esquerdo
    c_sibling = search(cie->get_cell_x(), cie->get_cell_y() + 1, cie->get_cell_level());
    if (c_sibling != NULL)
      L->push_back(c_sibling);
    
    //célula irmã do canto superior direito
    c_sibling = search(cie->get_cell_x() + 1, cie->get_cell_y() + 1, cie->get_cell_level());
    if (c_sibling != NULL)
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
  /*if (c1->get_cell_x() == c2->get_cell_x() &&
      c1->get_cell_y() == c2->get_cell_y() &&
      c1->get_cell_level() == c2->get_cell_level())*/
  if (c1 == c2)
    return true;
  return false;
}

list <cell *>::iterator mesh::merge (vector <cell *> * V, cell * c){
  list <cell *>::iterator ret;
  //bool allways_set_ret = false;

  //c->print_cell();
  
  //V->back() tem que ser do tipo "canto inferior esquerdo"
  char last_kind[2] = {((V->back())->get_cell_last_kind())[0], ((V->back())->get_cell_last_kind())[1]};
  assert(((V->back())->get_cell_kind())[0] == 'b' && ((V->back())->get_cell_kind())[1] == 'l');
  cell * new_c = new cell ((V->back())->get_cell_x() / 2, (V->back())->get_cell_y() / 2, last_kind, (V->back())->get_cell_level() - 1);
  if (new_c->get_cell_x() % 2 == 0 && new_c->get_cell_y() %2 == 0){
    last_kind[0] = 'b';
    last_kind[1] = 'l';
    new_c->set_cell_kind(last_kind);
  }
  else if (new_c->get_cell_x() % 2 == 1 && new_c->get_cell_y() % 2 == 0){
    last_kind[0] = 'b';
    last_kind[1] = 'r';
    new_c->set_cell_kind(last_kind);
  }
  else if (new_c->get_cell_x() % 2 == 0 && new_c->get_cell_y() % 2 == 1){
    last_kind[0] = 't';
    last_kind[1] = 'l';
    new_c->set_cell_kind(last_kind);
  }
  else{
    last_kind[0] = 't';
    last_kind[1] = 'r';
    new_c->set_cell_kind(last_kind);
  }
  //cout << "Nova célula: "; 
  //new_c->print_cell();
  insert(new_c);
  for (vector <cell *>::iterator it = V->begin(); it != V->end(); it++) {
    //(*it)->print_cell();
    if (!equal(*it, c)){
      remove(*it);
    }
  }
  ret = remove(c);
  //assert(allways_set_ret == true);
  return ret;
}
