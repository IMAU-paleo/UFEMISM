#! /bin/csh -f

cd src

make clean
make all

cd ..

rm -f UFEMISM_program

mv src/UFEMISM_program .
