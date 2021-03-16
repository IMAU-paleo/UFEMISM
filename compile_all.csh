#! /bin/csh -f

cd src_dev

#make clean
make all

cd ..

rm -f UFEMISM_program

mv src_dev/UFEMISM_program .