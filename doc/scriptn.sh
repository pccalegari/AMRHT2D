#Florianopolis fev 2024
# Nao esqueca de alterar o local de instalacao
# Update .bashrc file:
# Add: export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$path$/szip/lib:$path$/hdf5/lib:$path$/silo/lib:$path$/petsc/lib"

mkdir /home/priscila/Downloads/Programas
mkdir /home/priscila/Documents/Programas

MAIN_DOWNLOAD_DIR=/home/priscila/Downloads/Programas
MAIN_INSTALL_DIR=/home/priscila/Documents/Programas

wget -P $MAIN_DOWNLOAD_DIR https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz

wget -P $MAIN_DOWNLOAD_DIR https://zlib.net/zlib-1.3.1.tar.gz

wget -P $MAIN_DOWNLOAD_DIR https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar

wget -P $MAIN_DOWNLOAD_DIR https://github.com/LLNL/Silo/releases/download/4.11.1/silo-4.11.1.tar.xz

wget -P $MAIN_DOWNLOAD_DIR https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.20.4.tar.gz

cd $MAIN_DOWNLOAD_DIR

tar -xvzf szip-2.1.1.tar.gz
tar -xvzf zlib-1.3.1.tar.gz
tar -xvf hdf5-1.8.20.tar
tar -xvf silo-4.11.1.tar.xz
tar -xvzf petsc-3.20.4.tar.gz

cd szip-2.1.1
./configure --prefix=$MAIN_INSTALL_DIR/szip 
make 
make check
make install

cd ../zlib-1.3.1
./configure --prefix=$MAIN_INSTALL_DIR/zlib 
make 
make check
make install

cd ../hdf5-1.8.20
CC='cc -m64' ./configure --prefix=$MAIN_INSTALL_DIR/hdf5 --enable-fortran --enable-cxx --with-szlib=$MAIN_INSTALL_DIR/szip --with-zlib=$MAIN_INSTALL_DIR/zlib/include,$MAIN_INSTALL_DIR/zlib/lib
make 
make check
make install

cd ../silo-4.11.1
CC=gcc FC=gfortran CXX=g++ ./configure --prefix=$MAIN_INSTALL_DIR/silo --enable-fortran --with-szlib=$MAIN_INSTALL_DIR/szip --with-zlib=$MAIN_INSTALL_DIR/zlib/include/,$MAIN_INSTALL_DIR/zlib/lib/ --with-hdf5=$MAIN_INSTALL_DIR/hdf5/include/,$MAIN_INSTALL_DIR/hdf5/lib/
make 
make install

cd ../petsc-3.20.4
./configure --prefix=$MAIN_INSTALL_DIR/petsc --download-mpich --download-fblaslapack 
make PETSC_DIR=/home/priscila/Downloads/Programas/petsc-3.20.4 PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=/home/priscila/Downloads/Programas/petsc-3.20.4 PETSC_ARCH=arch-linux-c-debug install


rm -R MAIN_DOWNLOAD_DIR
