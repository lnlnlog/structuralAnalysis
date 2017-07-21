/*
 * Particle.cpp
 */

#include "Particle.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
using std::cout;
using std::endl;

Particle::Particle() { // constructor
    rx = ry = rz = 0.0; // zero xyz coordinates
    sigma = 0.0; // zero particle diameter
    id = -1; // assign id to be invalid
    clearVertex(); // clear vertex number
    isCrystalParticle = false; // assume this particle is liquid like with no info given
}

Particle::~Particle() { // destructor
    // empty
}

int Particle::setCrystal() { // label this particle as crystal like
    isCrystalParticle = true;
    return 0;
}

int Particle::unsetCrystal() { // label this particle as liquid like
    isCrystalParticle = false;
    return 0;
}


bool Particle::isCrystal() { // return the crystal or liquid label
    return isCrystalParticle;
}

int Particle::clearVertex() { // clear vp and vertex array
    vp = 0;
    for(int i=0; i<=MAXVERATOM - 1; i++)
        for(int j=0; j<=2; j++)
            vertex[i][j] = 0.0; // set vertex number to be invalid
    tx = ty = tz = 0.0;
    return 0;
}

int Particle::setPointT(double a, double b, double c) { // set the touching point
    tx = a;
    ty = b;
    tz = c;
    return 0;
}

int Particle::push_vertex(double x, double y, double z) { // add vertex coordinate to vertex array
    if(vp >= MAXVERATOM) {
        cout << "Error in Particle::push_vertex()" << endl;
        exit(1);
    }
    vertex[vp][0] = x;
    vertex[vp][1] = y;
    vertex[vp][2] = z;
    vp++;
    return 0;
}

bool Particle::isNeighbor() const { // true if this is a neighbor particle, vp > 0
    bool q;
    
    vp > 0 ? q = true : q = false;
    return q;
}

int Particle::print(ofstream& oFile, double cx, double cy, double cz) { // print vertices associated with this atom
    if(isNeighbor()) {
        sort(); // sort vertices so that they form clockwise order
        int i1;
        for(int i=0; i<=vp-1; i++) {
            i1 = (i + 1) % vp;
            oFile << vertex[i][0] + cx  << " " << vertex[i][1] + cy  << " " << vertex[i][2] + cz  << " "
                  << vertex[i1][0] + cx << " " << vertex[i1][1] + cy << " " << vertex[i1][2] + cz << endl;
        }
    }
    return 0;
}

double Particle::getarea() { // get voronoi polyhedron face area belongs to this neighbor
    double area = 0.0;
    
    if(isNeighbor()) {
        sort(); // sort vertices so that they form clockwise order
        int i1;
        for(int i=0; i<=vp-1; i++) {
            i1 = (i + 1) % vp;
            area += getTriangleArea(tx, ty, tz, vertex[i][0], vertex[i][1], vertex[i][2], vertex[i1][0], vertex[i1][1], vertex[i1][2]);
        }
    }
    else {
        area = 0.0;
    }
    
    return area;
}

double Particle::getvolumn() { // get pyramid volumn in order for the volumn of the entire voronoi polyhedron
    double area, h;
    
    area = getarea();
    h = sqrt(tx * tx + ty * ty + tz * tz);
    
    return area * h * 1.0 / 3.0;
}

double Particle::getTriangleArea(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2) const { // get area given coordinates of three vertices of triangle
    double a, b, c, p;
    
    a = sqrt(pow((x0 - x1),2.0) + pow((y0 - y1),2.0) + pow((z0 - z1),2.0));
    b = sqrt(pow((x1 - x2),2.0) + pow((y1 - y2),2.0) + pow((z1 - z2),2.0));
    c = sqrt(pow((x2 - x0),2.0) + pow((y2 - y0),2.0) + pow((z2 - z0),2.0));
    p = (a + b + c) / 2.0;
    
    double area;
    
    area = sqrt(p * (p - a) * (p - b) * (p - c));
    
    return area;
}

int Particle::sort() { //sort vertices into clockwise order
    int nv;
    double theta[MAXVERATOM];
    double a,b,c,r;
    double xi,yi,zi,ri,thetai;

    nv = vp;
    theta[0] = 0; // set angle of first vertex to be zero
    a = vertex[0][0] - tx;
    b = vertex[0][1] - ty;
    c = vertex[0][2] - tz;
    r = sqrt(a*a + b*b + c*c);
    for(int i=1; i<=nv-1; i++) { // search all the vertices but the first one
        xi = vertex[i][0] - tx;
        yi = vertex[i][1] - ty;
        zi = vertex[i][2] - tz;
        ri = sqrt(xi*xi + yi*yi + zi*zi);
        thetai = acos((a * xi + b * yi + c * zi) / r / ri);
        tx * (b * zi - c * yi) + ty * (c * xi - a * zi) + tz * (a * yi - b * xi) > 0 ? thetai *= 1.0 : thetai *= -1.0;
        theta[i] = thetai;
    }
    
    // sort according to theta[]
    int top = nv - 2;
    bool change = true;
    
    while(change && top >= 0) {
        change = false;
        
        for(int i=0; i<=top; i++) {
            int i1;
            i1 = i + 1;
            if(theta[i] > theta[i1]) {
                xi = vertex[i][0];
                yi = vertex[i][1];
                zi = vertex[i][2];
                thetai = theta[i];
                
                vertex[i][0] = vertex[i1][0];
                vertex[i][1] = vertex[i1][1];
                vertex[i][2] = vertex[i1][2];
                theta[i] = theta[i1];
                
                vertex[i1][0] = xi;
                vertex[i1][1] = yi;
                vertex[i1][2] = zi;
                theta[i1] = thetai;
                
                change = true;
            }
        }
        
        top--;
    }
    return 0;
}

int Particle::pop_vertex() { // pop one element from vertex array
    vp--;
    if(vp < 0) {
        cout << "Error in Particle::pop_vertex()" << endl;
        exit(1);
    }
    return 0;
}

int Particle::setrx(double x) { // set rx, return 0 if normal
    rx = x;
    return 0;
}

int Particle::setry(double y) { // set ry, return 0 if normal
    ry = y;
    return 0;
}

int Particle::setrz(double z) { // set rz, return 0 if normal
    rz = z;
    return 0;
}

int Particle::setsigma(double sig) { // set sigma, return 0 if normal
    sigma = sig;
    return 0;
}

int Particle::setid(int idnumber) { // set id, return 0 if normal
    id = idnumber;
    return 0;
}

