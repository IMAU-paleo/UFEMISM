#! /bin/csh -f

cd src_v1.1

make clean
make all

cd ..

rm -f UFEMISM_program

mv src_v1.1/UFEMISM_program .
