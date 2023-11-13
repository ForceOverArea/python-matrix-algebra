#!/bin/bash
echo [BUILD] building 'matrix.c'...
gcc -fPIC -c c_src/matrix.c -o c_src/matrix.o

echo [BUILD] building 'pymatlab_wrapper.c'...
gcc -fPIC -c c_src/pymatlab_wrapper.c -o c_src/pymatlab_wrapper.o

echo [BUILD] linking and building 'pymatlablib.so'...
gcc -fPIC -shared -o pymatlablib.so c_src/matrix.o c_src/pymatlab_wrapper.o

echo [BUILD] done!