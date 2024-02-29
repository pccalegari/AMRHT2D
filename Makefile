#ATENÇÃO:
# 1) Altere os caminhos onde foram feitas as instalaçoes 
# 2) NÃO ESQUEÇA DE ALTERAR O ARQUIVO .bashrc
#Inclua os caminhos onde foram feitas as instalaçoes 
#Adicione: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$path$/szip/lib:$path$/hdf5/lib:$path$/silo/lib:$path$/petsc/lib"


INCLUDE_PATH=include
LIBRARY_PATH=/home/priscila/Documentos/Programas
INCLUDE_SILO_PATH=$(LIBRARY_PATH)/silo/include
LIB_PATH=-L$(LIBRARY_PATH)/hdf5/lib -L$(LIBRARY_PATH)/szip/lib  -L$(LIBRARY_PATH)/zlib/lib -L$(LIBRARY_PATH)/silo/lib -lhdf5 -lsz -lz -lm -lsiloh5

PETSC_ARCH=""
PETSC_DIR=$(LIBRARY_PATH)/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_FLAGS = -I${PETSC_DIR}/include -L${PETSC_DIR}/lib/libpetsc_real.so

SRC=src
TEST=test
DATA=data
CC=g++

OPTIONS=-Wall -ansi -pedantic -Wno-unused-result -O3 -std=c++11
LIBS= -L$(LIBRARY_PATH)/silo/lib -lsiloh5 -L$(LIBRARY_PATH)/hdf5/lib -lhdf5 -L$(LIBRARY_PATH)/zlib/lib -lz -lm

all:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/point.cpp -lm -o $(SRC)/point.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o
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
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o	
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test.cpp -lm -o $(TEST)/test.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/weight.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test1

test2:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_particles.cpp -lm -o $(TEST)/test.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/weight.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test2

test3:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c $(PETSC_FLAGS) -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_Laplace.cpp -lm -o $(TEST)/test.o ${PETSC_KSP_LIB}
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/weight.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test3 ${PETSC_VEC_LIB}

test4:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/point.cpp -lm -o $(SRC)/point.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/mls.cpp -lm -o $(SRC)/mls.o

	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c $(PETSC_FLAGS) -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_LaplaceAMR.cpp -lm -o $(TEST)/test.o ${PETSC_KSP_LIB}
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/weight.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test4 ${PETSC_VEC_LIB}

test5:
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/particle.cpp -lm -o $(SRC)/particle.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/ode.cpp -lm -o $(SRC)/ode.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/weight.cpp -lm -o $(SRC)/weight.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c $(PETSC_FLAGS) -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/exemplo_QR.cpp -lm -o $(TEST)/test.o ${PETSC_KSP_LIB}
	$(CC) -fPIC $(OPTIONS)  $(SRC)/particle.o $(SRC)/ode.o $(SRC)/weight.o $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test.o $(LIBS) -o $(TEST)/test5 ${PETSC_VEC_LIB}


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

runtest3:
	./test/test3	

runtest4:
	./test/test4

cleansilo:
	rm -rf $(DATA)/*.silo
