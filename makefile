ising: ising.o Lattice.o
	g++ ising.o Lattice.o -o ising
ising.o: ising.cpp ising.hpp
	g++ -c ising.cpp -o ising.o
Lattice.o: Lattice.cpp Lattice.hpp
	g++ -c Lattice.cpp -o Lattice.o
clean:
	rm -r *.o
clear:
	rm -r *.dat IMAGES/*.pgm
