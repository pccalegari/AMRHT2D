#ATENÇÃO: NÃO ESQUEÇA DE ALTERAR O ARQUIVO .bashrc
#Adicione: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/szip/lib:/opt/hdf5/lib:/opt/silo/lib"
CC=g++
#CC=/opt/hdf5/bin/h5c++
INCLUDE_PATH=../include 
INCLUDE_SILO_PATH=/opt/silo/include/
#/opt/silo/include

#INCLUDE_HDF5_PATH=/opt/hdf5/include
LIB_PATH=-L/opt/hdf5/lib -L/opt/szip/lib -L/opt/zlib/lib -L/opt/silo/lib -lhdf5 -lsz -lz -lm -lsiloh5
#-L/opt/hdf5/lib -L/opt/szip/lib  -L/opt/zlib/lib -L/opt/silo/lib -lhdf5 -lsz -lz -lm -lsiloh5
SRC=../src
TEST=../test
#-Wl,--rpath -Wl,LIBDIR

OPTIONS=-Wall -ansi -pedantic -Wno-unused-result -O3 -std=c++11
#LIBS= /opt/hdf5/lib -lhdf5 -L/opt/szip/lib -lsz  -L/opt/zlib/lib -lz -L/opt/silo/lib -lsiloh5 -lm
#LIBS=-L/opt/zlib/lib/libz.so -L/opt/szip/lib/libsz.so -L/opt/hdf5/lib/libhdf5.so -L/opt/silo/lib/libsiloh5.la -lm
LIBS=-L/opt/silo/lib -lsiloh5 -L/opt/hdf5/lib -lhdf5 -L/opt/zlib/lib -lz -lm
#-L/opt/silo/lib -lsiloh5 -L/opt/hdf5/lib -lhdf5 -L/opt/zlib/lib -lz -lm
all:  
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/cell.cpp -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) $(SRC)/dominio.cpp -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test_triangular_mesh.cpp -lm -o $(TEST)/test_triangular_mesh.o
	$(CC) -fPIC $(OPTIONS)  $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test_triangular_mesh.o $(LIBS) -o $(TEST)/test_triangular_mesh

gdb:
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/cell.cpp $(LIBS) -lm -o $(SRC)/cell.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/double_cell.cpp $(LIBS) -lm -o $(SRC)/double_cell.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/hash_table.cpp $(LIBS) -lm -o $(SRC)/hash_table.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/double_hash_table.cpp $(LIBS) -lm -o $(SRC)/double_hash_table.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) $(SRC)/dominio.cpp $(LIBS) -lm -o $(SRC)/dominio.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(SRC)/mesh.cpp $(LIBS) -lm -o $(SRC)/mesh.o
	$(CC) $(OPTIONS) -c -g -I$(INCLUDE_PATH) -I$(INCLUDE_SILO_PATH) $(TEST)/test.cpp $(LIBS) -lm -o $(TEST)/test_triangular_mesh.o
	$(CC) $(OPTIONS) -g $(SRC)/cell.o $(SRC)/double_cell.o $(SRC)/hash_table.o $(SRC)/double_hash_table.o $(SRC)/dominio.o $(SRC)/mesh.o $(TEST)/test_triangular_mesh.o $(LIB_PATH) $(LIBS) -lm -o $(TEST)/test_triangular_mesh	

clean:
	rm -rf $(SRC)/*.o $(TEST)/*.o $(TEST)/test_triangular_mesh   


