gfortran  -c -w -fPIC -fopenmp Wrapper.f90
gcc -c -w -fPIC -fopenmp Parallel.c 
gcc -shared -fopenmp -o libSofia.so Parallel.o Wrapper.o
