/*
 * Particle.hpp
 * declare Particle class in voronoi tessellation
 *
 */

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <fstream>

using std::ofstream;
using std::ifstream;

const int MAXVERATOM = 20; // max number of vertices one atom could have, depending on dimension, MAXVERATOM = 2 in 2 dimension

class Particle {
private:
    double rx,ry,rz; // xyz coordinates
    double sigma; // atom diameter
    int id; // particle id 1..n species
    int vp; // vertex pointer
    double vertex[MAXVERATOM][3]; // vertex number associate with this particle
    double tx,ty,tz; // point T coordinate, joint of a plane and perpendicular line(starts from origin)
    bool isCrystalParticle; // true if this particle is crystal like
    
    int sort(); //sort vertices into clockwise order
    double getTriangleArea(double, double, double, double, double, double, double, double, double) const; // get area given coordinates of three vertices of triangle


public:
    Particle(); // constructor
    ~Particle(); // destructor
    double getrx() const { // get rx
        return rx;
    }
    double getry() const { // get ry
        return ry;
    }
    double getrz() const { // get rz
        return rz;
    }
    double getsigma() const { // get sigma
        return sigma;
    }
    int getid() const { // get id
        return id;
    }
    int getvp() const { // get vp
        return vp;
    }
    int setrx(double ); // set rx, return 0 if normal
    int setry(double ); // set ry, return 0 if normal
    int setrz(double ); // set rz, return 0 if normal
    int setsigma(double ); // set sigma, return 0 if normal
    int setid(int ); // set id, return 0 if normal
    int clearVertex(); // clear vp and vertex array
    int push_vertex(double, double, double); // add vertex coordinate to vertex array
    int pop_vertex(); // pop one vertex number from vertex array
    int print(ofstream&, double, double ,double); // print vertices associated with this atom
    bool isNeighbor() const; // true if this is a neighbor particle, vp > 0
    int setPointT(double, double, double); // set the touching point
    double getarea(); // get voronoi polyhedron face area belongs to this neighbor
    double getvolumn(); // get pyramid volumn in order for the volumn of the entire voronoi polyhedron
    int setCrystal(); // label this particle as crystal like
    int unsetCrystal(); // label this particle as liquid like
    bool isCrystal(); // return the crystal or liquid label


};

#endif
