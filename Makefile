all: main.o Particle.o Box.o Q6.o
	icpc -o voronoi main.o Particle.o Box.o Q6.o
main.o:
	icpc -g -O3 -Wall -c main.cpp
Particle.o:
	icpc -g -O3 -Wall -c Particle.cpp
Box.o:
	icpc -g -O3 -Wall -c Box.cpp
Q6.o:
	icpc -g -O3 -Wall -c Q6.cpp
cleanup:
	rm *.o
	rm voronoi
