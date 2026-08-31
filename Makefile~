#ATENÇÃO: NÃO ESQUEÇA DE ALTERAR O ARQUIVO .bashrc
#Adicione: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/szip/lib:/opt/hdf5/lib:/opt/silo/lib"
CC=g++

INCLUDE_PATH=include 
INCLUDE_SILO_PATH=/opt/silo/include
LIB_PATH=-L/opt/hdf5/lib -L/opt/szip/lib  -L/opt/zlib/lib -L/opt/silo/lib -lhdf5 -lsz -lz -lm -lsiloh5
SRC=src
TEST=test
DATA=data

OPTIONS=-Wall -ansi -pedantic -Wno-unused-result -O3 -std=c++11
LIBS= -L/opt/silo/lib -lsiloh5 -L/opt/hdf5/lib -lhdf5 -L/opt/zlib/lib -lz -lm

all:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test.cpp -lm -o $(TEST)/test.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test

test1:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test.cpp -lm -o $(TEST)/test.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test1

test2:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_particles.cpp -lm -o $(TEST)/test.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test2

gdb:
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/cell.cpp $(LIBS) -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp $(LIBS) -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp $(LIBS) -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp $(LIBS) -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/dominio.cpp $(LIBS) -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp $(LIBS) -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_particles.cpp $(LIBS) -lm -o $(TEST)/test.o
	$(CC) $(OPTIONS) -g $(SRC)/particle.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIB_PATH) $(LIBS) -lm -o $(TEST)/test	

clean:
	rm -rf $(SRC)/*.o $(TEST)/*.o $(TEST)/test $(DATA)/*.silo  

run:
	./test/test

runtest1:
	./test/test1

runtest2:
	./test/test2

cleansilo:
	rm -rf $(DATA)/*.silo
