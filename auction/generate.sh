mpiicc -c floauction.c -openmp
mpiifort -g -O3 -cpp -o generateRandMat generateRandMat.f90 floauction.o -openmp 
