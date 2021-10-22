export FC=gfortran
cmake -DSERIAL=ON ../src
cmake --build . --parallel
