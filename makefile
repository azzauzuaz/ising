ising: ising.o Lattice.o Histogram.o
	g++ ising.o Lattice.o Histogram.o -o ising
ising.o: ising.cpp ising.hpp
	g++ -c ising.cpp -o ising.o
Lattice.o: Lattice.cpp Lattice.hpp
	g++ -c Lattice.cpp -o Lattice.o
Histogram.o: Histogram.cpp Histogram.hpp
	g++ -c Histogram.cpp -o Histogram.o
clean:
	rm -r *.o
