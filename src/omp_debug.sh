gfortran -g -c -w -fPIC Wrapper.f90
gcc -g -c -w -fPIC Parallel.c 
gcc -g -shared -o libSofia.so Parallel.o Wrapper.o -lgfortran
