/* 
 * Box.hpp
 * declare class Box
 */

#ifndef BOX_HPP
#define BOX_HPP

#include "Particle.hpp"
#include <complex>

using std::complex;

const int MaxAtomNumber = 9000; // max of npart allowed
const int MAXCAN = 500; // max of candidates in voron
const int MAXVER = 500; // max of vertices of a given atom
const double TOL = 1e-10; // tolerance in det

class Box {
private:
    Particle atom[MaxAtomNumber]; // atom info
    double Rratio; // bidisperse ratio
    double canShellInput; // general shell radius, need to be specified with box size
    double canShell; // shell radius used in each frame
    int npart; // number of particles
    double boxx,boxy,boxz; // box dimension
    double boxx_last;
    int frameStart, frameInter; // frame controller
    char movieFile[50]; // file from which we extract configs
    char boxFile[50]; // file from which we extract box length
    double composition; // small particle number fraction
    complex<double> qlm[MaxAtomNumber][13]; // for q6 only, 13 = 2 * 6 + 1
    
    int initializeFromFile(ifstream&, ifstream&, bool&); // read from file to initialize atom[]
    int printFrame() const; // print frames for visualization
    double inBoxx(double ) const; // pbc in x direction
    double inBoxy(double ) const; // pbc in y direction
    double inBoxz(double ) const; // pbc in z direction
    int sort(int, double*, double*, double*, double*, int*, int&); // sort neighbors into increasing distance order
    int work(int, int, double*, double*, double*, double*, int*, bool&); // function performing voronoi tessellation
    int voron(int); // main loop doing voronoi tessellation    
    int surrounding(int [MaxAtomNumber][MaxAtomNumber], complex<double> [MaxAtomNumber][13], int, double*, bool*);
    
    int makeUnion(int, int, int [], int[]); //union operation in dynamic connectivity
    int root(int ,int []);
    double distance(int, int) const; // distance between two particles
    int min2(double&, double&); // min function
    int max2(double&, double&); // max function
    
public:
    Box(); // constructor
    ~Box(); // destructor
    int run(); // read couple of frames and voronoi tessellation
};

#endif
