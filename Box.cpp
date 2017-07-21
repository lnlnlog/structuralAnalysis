/* 
 * Box.cpp
 */

#include "Box.hpp"
#include "const.hpp"
#include "Q6.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <string>

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setprecision;
using std::string;
using std::streampos;

Box::Box() { // constructor
    
    cout << "Enter npart:" << endl;
    cin  >> npart;
    cout << "Enter Ratio:" << endl;
    cin  >> Rratio;
    cout << "Enter canShellInput:" << endl;
    cin  >> canShellInput;
    cout << "Enter movie file name:" << endl;
    cin  >> movieFile;
    cout << "Enter box file name:" << endl;
    cin  >> boxFile;
    cout << "Enter frame number to stop:" << endl;
    cin  >> frameStart;
    cout << "Enter frame interval:" << endl;
    cin  >> frameInter;
    cout << "Enter composition:" << endl;
    cin  >> composition;
    
    boxx = boxy = boxz = 0.0; // initialize box dimension

    for(int i=0;i<=MaxAtomNumber-1;i++)
        for(int j=0;j<=13;j++)
            qlm[i][j] = complex<double> (0.0,0.0);
}

Box::~Box() { // destructor
    // empty
}

int Box::run() { // read couple of frames and voronoi tessellation
    
    ifstream fileMovie,fileBox;
    
    fileMovie.open(movieFile, ios::in);
    fileBox.open(boxFile, ios::in);
    
//    double a,b,c;
//    streampos begins,ends;
//    
//    begins = fileBox.tellg();
//    fileBox >> a >> b >> c;
//    ends   = fileBox.tellg();
//    
//    fileBox.seekg(-(ends-begins),ios::end);
//    fileBox >> a >> b >> c;
//    boxx_last = a;
//    
//    cout << "boxx_last = " << boxx_last << endl; // may be problematic
//        
//    fileBox.seekg(0,ios::beg);
    
    bool fileEnd = false;
    char ch;
    double x, y, z;
    
    for(int frame = 1; frame < frameStart; frame++) {
        fileMovie >> npart;
        for(int i=0; i < npart; i++)
            fileMovie >> ch >> x >> y >> z;
        
        fileBox >> x >> y >> z;
    }
        
    
    for(int frame = frameStart; ; frame++) {
        
        initializeFromFile(fileMovie, fileBox, fileEnd);
        
        if(fileEnd) break;
    
        if( (frame - frameStart) % frameInter == 0 ) {
            
            double sigma;
            sigma = 1.0;
            
            canShell = canShellInput * sigma * Rratio;
            
            voron(frame);
            
            cout << "frame = " << frame << endl;
            cout << std::setprecision(10) << "boxx = " << boxx << endl;
            cout << std::setprecision(10) << "boxy = " << boxy << endl;
            cout << std::setprecision(10) << "boxz = " << boxz << endl;
        }

    }
    
    cout << std::setprecision(10) << "boxx = " << boxx << endl;
    cout << std::setprecision(10) << "boxy = " << boxy << endl;
    cout << std::setprecision(10) << "boxz = " << boxz << endl;
    
    fileMovie.close();
    fileBox.close();
    
    return 0;
}

int Box::initializeFromFile(ifstream& fileIn, ifstream& fileInBox, bool& fileEndReached) { // read from file to initialize atom[]
    char ch;
    double x,y,z;
    
    fileIn >> npart;
    
    if(npart > MaxAtomNumber) {
        cout << "MaxAtomNumber reached" << endl;
        exit(0);
    }
    
    if(fileIn.eof()) {
        fileEndReached = true;
        return 0;
    }
    
    int nsLimit;
    
    nsLimit = static_cast<int> (round(npart * composition));
    
//    boxx = boxy = boxz = pow(PI * (nsLimit + pow(Rratio,3.0) * (npart - nsLimit)) / 6.0 / phi,1.0/3.0);
    
    fileInBox >> boxx >> boxy >> boxz;
    
    for(int i = 0; i<=npart - 1; i++) {
        fileIn >> ch >> x >> y >> z;
        
        x = inBoxx(x); // pbc on
        y = inBoxy(y); // pbc on
        z = inBoxz(z); // pbc on

        if(ch == 'N') {
            atom[i].setid(1);
            atom[i].setsigma(1.0); // assume in movie file, atom diameter is 1.0
        }
        else if(ch == 'O') {
            atom[i].setid(2);
            atom[i].setsigma(1.0 * Rratio); // large particle diameter
        }
        else {
            cout << "error initializeFromFile: wrong atom character" << endl;
            exit(1);
        }
        atom[i].setrx(x);
        atom[i].setry(y);
        atom[i].setrz(z);
    }
    
    return 0;
}

int Box::printFrame() const { // print frames for visualization
    ofstream outFile;
    char ch;
    outFile.open("config.txt",ios::out);
    
    outFile << npart << endl;
    outFile << boxx << " " << boxy << " " << boxz << endl;
    for(int i=0; i<= npart-1; i++) {
        atom[i].getid() == 1 ? ch = 'N' : ch = 'O';
        outFile << ch << " " << inBoxx(atom[i].getrx()) << " " << inBoxy(atom[i].getry()) << " " << inBoxz(atom[i].getrz()) << endl;
    }
    
    outFile.close();
    return 0;
}

int Box::voron(int framenum) { // main loop doing voronoi tessellation
    double rxj,ryj,rzj;
    int can,ncan;
    double rxij,ryij,rzij,rijsq;
    double canShellsq;
    double px[MAXCAN],py[MAXCAN],pz[MAXCAN],ps[MAXCAN];
    int tag[MAXCAN];
    int QQn[MaxAtomNumber][MaxAtomNumber];
    
    bool isOpen[MaxAtomNumber];
    double bop[MaxAtomNumber];
    
    ofstream oFile,oFileConfig1,oFileConfig2,oFileQ6;
//    oFile.open("packingFNumber.txt",ios::out);
//    oFileConfig1.open("shell.txt",ios::out);
//    oFileConfig2.open("vertexConfig2.txt",ios::app);
//    oFileQ6.open("surfBulk.txt",ios::app);
    
    canShellsq = canShell * canShell;
    
    double volumeTotal = 0.0;
    double particleVolumeTotal = 0.0;
    
    for(int j = 0; j <= npart-1; j++) {
        rxj = atom[j].getrx();
        ryj = atom[j].getry();
        rzj = atom[j].getrz();
        can = 0; // px,py,pz,ps,tag start from 0
        for(int i = 0; i <= npart - 1; i++) {
            if(i == j) continue;
            rxij = atom[i].getrx() - rxj;
            ryij = atom[i].getry() - ryj;
            rzij = atom[i].getrz() - rzj;
            
            rxij = inBoxx(rxij);
            ryij = inBoxy(ryij);
            rzij = inBoxz(rzij);
            
            rijsq = rxij * rxij + ryij * ryij + rzij * rzij;
            
            if(rijsq < canShellsq) {
                can++;
                if(can > MAXCAN) {
                    cout << "Error in Box::voron(): too many candidates" << endl;
                    cout << "can = " << can << "MAXCAN = " << MAXCAN << endl;
                    exit(1);
                }
                px[can - 1] = rxij;
                py[can - 1] = ryij;
                pz[can - 1] = rzij;
                ps[can - 1] = rijsq;
                tag[can - 1] = i;
            }
        } //i, candidates selected
        
        ncan = can; // number of candidates
        
        sort(j, px, py, pz, ps, tag, ncan);

        work(j, ncan, px, py, pz, ps, tag, isOpen[j]);
        
        double wx[MAXCAN], wy[MAXCAN], wz[MAXCAN]; // bond vector info with respect to central atom
        double areaCan[MAXCAN]; // voronoi polyhedron, face area belongs to one neighbor
        int nnb = 0; // number of nearest neighbor
        double volumeJ = 0.0; // volume for central particle
        
//        if(isOpen[j] == false) { // if bulk particle
            for(int k=0; k <= ncan-1; k++) {
                if(atom[tag[k]].isNeighbor()) { // print out one face of polygon
                    wx[nnb] = px[k];
                    wy[nnb] = py[k];
                    wz[nnb] = pz[k];
//                    areaCan[nnb] = 1;//atom[tag[k]].getarea(); // get voronoi polyhedron face area belongs to this neighbor
                    if(isOpen[j] == false) areaCan[nnb] = atom[tag[k]].getarea(); // for a bulk particle
                    else areaCan[nnb] = 1; // for a surface particle
                    
                    if(isOpen[j] == false) volumeJ += atom[tag[k]].getvolumn();
                    
                    QQn[j][nnb] = tag[k];
//
                    nnb++;
//                    atom[tag[k]].print(oFile, rxj, ryj, rzj); // print out all edges
                }
            }
        
            if(isOpen[j] == false) {
                volumeTotal += volumeJ;
                particleVolumeTotal += PI / 6.0 * pow(atom[j].getsigma(),3.0);
//                oFile << PI / 6.0 * pow(atom[j].getsigma(),3.0) / volumeJ << endl;
            }
        
            QQn[j][nnb] = -1;
        
            Q6::w(6, wx, wy, wz, areaCan, nnb, qlm, j);
            
            // Q6, qlm[j][0..12]
//            bop[j] = 0.0;
//            for(int i=0;i<=12;i++) bop[j] += norm(qlm[j][i]);
//            
//            bop[j] = sqrt(bop[j] * 4.0 * PI / (2 * 6 + 1));
        
//            oFileQ6 << bop[j] << endl;
//        }
//        else { // if surface particle
//            QQn[j][nnb] = -1;
//            for(int i=0;i<=12;i++) qlm[j][i] = complex<double> (0.0,0.0);
//            bop[j] = 0.0;
//        }
        
//        oFileQ6 << bop[j] << endl;
        
    } // j
    
    char ch;
    
//    oFileConfig1 << npart << endl << endl;
//    oFileConfig2 << npart << endl << endl;
    
    int ac_S = 0, ac_B = 0;
    for(int i=0;i<=npart-1;i++) {
        if(isOpen[i]) ch = 'O';
        else ch = 'N';
        
        if(isOpen[i]) {
            ac_S++;
//            oFileConfig1 << ch << " " << atom[i].getrx() << " " << atom[i].getry() << " " << atom[i].getrz() << endl;
//            oFileConfig2 << 'N' << " " << 100 << " " << 0 << " " << 0 << endl;
        }
        else {
            ac_B++;
//            oFileConfig1 << ch << " " << atom[i].getrx() << " " << atom[i].getry() << " " << atom[i].getrz() << endl;
//            oFileConfig2 << 'N' << " " << atom[i].getrx() << " " << atom[i].getry() << " " << atom[i].getrz() << endl;
        }
    }
    
    cout << "surface particle = " << ac_S << " " << "Bulk particle = " << ac_B << endl;
    
//    cout << "packing fraction = " << particleVolumeTotal / volumeTotal << endl;
//    cout << "number density   = " << ac_B / volumeTotal << endl;
    
//    oFile << particleVolumeTotal / volumeTotal << " " << ac_B / volumeTotal << endl;
    
//    oFileQ6 << framenum << " " << ac_S << " " << ac_B << endl;
    
    surrounding(QQn,qlm,framenum,bop,isOpen);
//    oFile.close();
//    oFileConfig1.close();
//    oFileConfig2.close();
//    oFileQ6.close();
    return 0;
}

int Box::surrounding(int Qn[MaxAtomNumber][MaxAtomNumber], complex<double> qm[MaxAtomNumber][13], int fnumber, double* bopNumber, bool* openPart) {
    
    string ch;
    
    const int l = 6; // Q6, l = 6
    int nOfCrystalPart = 0;
    
    ofstream q1q2,qConfigB,qConfigS,qConfig,qC;
    q1q2.open("nOfCrystal.txt",ios::app);
//    qConfig.open("cluster.txt",ios::out);
//    qC.open("cluster4.txt",ios::out);
//    qConfigB.open("crystalConfigB.txt",ios::app);
//    qConfigS.open("crystalConfigS.txt",ios::app);
    
//    qConfigB << npart << endl << endl;
//    qConfigS << npart << endl << endl;
    
    for(int j=0; j<npart; j++) {
        int i=0;
        int nconn = 0;
        atom[j].unsetCrystal(); // assume this is liquid like particle with no info given
        
        if(openPart[j] == false) { // if bulk particle
            while(Qn[j][i] != -1) {
                double sum = 0.0, d1 = 0.0, d2 = 0.0;
                
//                if(openPart[Qn[j][i]] == false) {
                    for(int m=-l; m<=l; m++) {
                        sum += (qm[j][m + l] * conj(qm[Qn[j][i]][m + l])).real();
                        d1 += norm(qm[j][m + l]);
                        d2 += norm(qm[Qn[j][i]][m + l]);
                    }
                    
                    if(d1 > 0.0 && d2 > 0.0) sum = sum / pow((d1 * d2),0.5);
//                }
                
                if(sum > 0.7) nconn++;
                
                i++;
            }
        
//            q1q2 << nconn << endl;
        
            if(nconn >= 7) {
                ch = "N";
//                qConfigB << ch << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
                nOfCrystalPart++;
                atom[j].setCrystal();
            }
            else {
//                ch = "N";
//                qConfigB << ch << " " << 100 << " " << 0 << " " << 0 << endl;
                  ch = "O";
//                  qConfigB << ch << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
            }
//            qConfigS << "O" << " " << 100 << " " << 0 << " " << 0 << endl; // take a position so S,B files same length
        }
//        else { // if surface particle
//            while(Qn[j][i] != -1) {
//                double sum = 0.0, d1 = 0.0, d2 = 0.0;
//                
////                if(openPart[Qn[j][i]] == true) {
//                    for(int m=-l; m<=l; m++) {
//                        sum += (qm[j][m + l] * conj(qm[Qn[j][i]][m + l])).real();
//                        d1 += norm(qm[j][m + l]);
//                        d2 += norm(qm[Qn[j][i]][m + l]);
//                    }
//                    
//                    if(d1 > 0.0 && d2 > 0.0) sum = sum / pow((d1 * d2),0.5);
////                }
//                
//                if(sum > 0.7) nconn++;
//                
//                i++;
//            }
//            
////            q1q2 << nconn << endl;
//            
//            if(nconn >= 5) {
////                ch = "O";
////                qConfigS << ch << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
//                ch = "N";
////                qConfigS << ch << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
//
//                nOfCrystalPart++;
//                atom[j].setCrystal();
//            }
//            else {
////                ch = "O";
////                qConfigS << ch << " " << 100 << " " << 0 << " " << 0 << endl;
//                ch = "O";
////                qConfigS << ch << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
//            }
////            qConfigB << "N" << " " << 100 << " " << 0 << " " << 0 << endl; // take a position so S,B files same length
//        }
        
        
    }
    
    cout << "nOfCrystalPart = " << nOfCrystalPart << endl;
    
    int nOfCrystalPart_surf = 0;
    
    for(int j=0;j<npart;j++) {
//    for(int j=4;j<=4;j++) {

        if(openPart[j] == true) { // if a surface particle
            int i = 0;
            int bulkNeighborCount = 0;
            
            while(Qn[j][i] != -1) { // loop over its neighbors
//                if(atom[Qn[j][i]].isCrystal() && openPart[Qn[j][i]] == false) { // if a bulk crystal particle
                    bulkNeighborCount++;
//                }
                i++;
            }
            
//            cout << "bulkNeighborCount = " << bulkNeighborCount << endl;
            
            if(bulkNeighborCount <= 2) {
                atom[j].unsetCrystal();
            }
            else { // bulkNeighborCount >= 3
                
                double shapeNumber = 1000; // initialize it with a large number
                
                for(int i1 = 0;i1 <= i-3;i1++) {
//                    if(atom[Qn[j][i1]].isCrystal() && openPart[Qn[j][i1]] == false) {
                        for(int i2 = i1 + 1;i2 <= i-2;i2++) {
//                            if(atom[Qn[j][i2]].isCrystal() && openPart[Qn[j][i2]] == false) {
                                for(int i3 = i2 + 1;i3 <= i-1;i3++) {
//                                    if(atom[Qn[j][i3]].isCrystal() && openPart[Qn[j][i3]] == false) {
//                                        q1q2 << distance(j,Qn[j][i1]) << " "
//                                             << distance(j,Qn[j][i2]) << " "
//                                             << distance(j,Qn[j][i3]) << " "
//                                             << distance(Qn[j][i1],Qn[j][i2]) << " "
//                                             << distance(Qn[j][i1],Qn[j][i3]) << " "
//                                             << distance(Qn[j][i2],Qn[j][i3]) << endl;
//                                        cout << distance(j,Qn[j][i1]) << " "
//                                        << distance(j,Qn[j][i2]) << " "
//                                        << distance(j,Qn[j][i3]) << " "
//                                        << distance(Qn[j][i1],Qn[j][i2]) << " "
//                                        << distance(Qn[j][i1],Qn[j][i3]) << " "
//                                        << distance(Qn[j][i2],Qn[j][i3]) << endl;
//                                        cout << i1 << " " << i2 << " " << i3 << endl;
                                    
                                        double dis[6],mindis,maxdis;
                                        
                                        dis[0] = distance(j,Qn[j][i1]);
                                        dis[1] = distance(j,Qn[j][i2]);
                                        dis[2] = distance(j,Qn[j][i3]);
                                        dis[3] = distance(Qn[j][i1],Qn[j][i2]);
                                        dis[4] = distance(Qn[j][i1],Qn[j][i3]);
                                        dis[5] = distance(Qn[j][i2],Qn[j][i3]);
                                        
                                        mindis = dis[0];
                                        maxdis = dis[0];
                                        
                                        for(int i4 = 1;i4 <= 5;i4++) {
                                            min2(mindis,dis[i4]);
                                            max2(maxdis,dis[i4]);
                                        }

                                        maxdis = maxdis - mindis;
                                        min2(shapeNumber,maxdis);
                                    
//                                    }
                                }
//                            }
                        }
//                    }
                }
                    
//                q1q2 << j << " " << shapeNumber << endl;
                
                if(shapeNumber < 0.0482) nOfCrystalPart_surf++;
            } // endif bulkneighborcount
            
//            qC << 3 << endl << endl;
//            
//            qC << 'N' << " " << atom[Qn[j][0]].getrx() << " " << atom[Qn[j][0]].getry() << " " << atom[Qn[j][0]].getrz() << endl;
//            qC << 'N' << " " << atom[Qn[j][2]].getrx() << " " << atom[Qn[j][2]].getry() << " " << atom[Qn[j][2]].getrz() << endl;
//            qC << 'N' << " " << atom[Qn[j][4]].getrx() << " " << atom[Qn[j][4]].getry() << " " << atom[Qn[j][4]].getrz() << endl;

            
//                qConfig << bulkNeighborCount+1 << endl << endl;
//            
//                i = 0;
//            
//                qConfig << 'O' << " " << atom[j].getrx() << " " << atom[j].getry() << " " << atom[j].getrz() << endl;
//                while(Qn[j][i] != -1) { // loop over its neighbors
//                    cout << "crystal " << atom[Qn[j][i]].isCrystal() << " surface " << openPart[Qn[j][i]] << endl;
//            
////                    if(atom[Qn[j][i]].isCrystal() && openPart[Qn[j][i]] == false)
//                    qConfig << 'N' << " " << atom[Qn[j][i]].getrx() << " " << atom[Qn[j][i]].getry() << " " << atom[Qn[j][i]].getrz() << endl;
//                    i++;
//                }
            
        }
    }

    cout << "nOfCrystalPart = " << nOfCrystalPart_surf << endl;
    q1q2 << fnumber << " " << nOfCrystalPart_surf << " " << nOfCrystalPart << endl;
    
    // here begins calculation of number of clusters
    
//    int id[MaxAtomNumber];
//    int sz[MaxAtomNumber];
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            id[j] = j;
//            sz[j] = 1;
//        }
//        else {
//            id[j] = -1;
//            sz[j] = -1;
//        }
//    }
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            int i=0;
//            while(Qn[j][i] != -1) {
//                if(atom[Qn[j][i]].isCrystal()) {
//                    makeUnion(j,Qn[j][i],id,sz);
//                }
//                i++;
//            }
//        }
//    }
    
//    int clusterId[MaxAtomNumber] = {0};
//    
//    int id2[MaxAtomNumber];
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            id2[j] = root(j,id);
//        }
//    }
//    for(int j=0; j<npart; j++) {
//        id[j] = id2[j];
//    }
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            clusterId[id[j]]++;
//        }
//    }
//    
//    int nCluster = 0;
//    for(int j=0; j<npart; j++) {
//        if(clusterId[j] >= 13) { // define cluster size > 13 as a cluster, 13 default
//            nCluster++;
//        }
//        else if(clusterId[j] > 0){
//            for(int k=0; k<npart; k++) {
//                if(id[k] == j) atom[k].unsetCrystal();
//            }
//            nOfCrystalPart -= clusterId[j];
//        }
//    }
    
//    q1q2 << fnumber << " " << static_cast<double>(nOfCrystalPart) / static_cast<double>(npart) << endl;
    
    //{"Li","C","N","O","F","Br","I","He","P","S","B","H","Be","Ti","Fe"}
//    string cheElement[15] = {"Li","Ti","N","O","F","Br","I","He","P","S","B","H","Be","Fe","C"};
//    int intElement[15] = {-1};
//    int len = 0;
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            int i=0;
//            for(i=0; i<len; i++) {
//                if(id[j] == intElement[i]) break;
//            }
//            if(i == len) {
//                intElement[len] = id[j];
//                len++;
//            }
//        }
//    }
    
    // cluster size
//    for(int j=0; j<len; j++) {
//        cout << "size[" << j << "] = " << sz[intElement[j]] << endl;
//    }
//    
//    for(int j=0; j<npart; j++) {
//        if(atom[j].isCrystal()) {
//            int i=0;
//            for(i=0; i<len; i++) {
//                if(id[j] == intElement[i]) break;
//            }
//            if(i==len) {
//                cout << "Error in print" << endl;
//                exit(0);
//            }
    
//            q1q2 << cheElement[i] << " " << inBoxx(atom[j].getrx()) << " " << inBoxy(atom[j].getry()) << " " << inBoxz(atom[j].getrz()) << endl;
//        }
//    }
    
    q1q2.close();
//    qConfig.close();
//    qC.close();
//    qConfigB.close();
//    qConfigS.close();
    
    return 0;
}

int Box::min2(double& a, double& b) {
    
    if(b < a) a = b;
    
    return 0;
}

int Box::max2(double& a, double& b) {
    
    if(b > a) a = b;
    
    return 0;
}

double Box::distance(int a, int b) const { // distance between two particles
    double rxij,ryij,rzij,r,sigma;
    
    sigma = (atom[a].getsigma() + atom[b].getsigma()) / 2.0;
    
    rxij = atom[a].getrx() - atom[b].getrx();
    ryij = atom[a].getry() - atom[b].getry();
    rzij = atom[a].getrz() - atom[b].getrz();
    
    rxij = inBoxx(rxij);

    r = sqrt(rxij * rxij + ryij * ryij + rzij * rzij);

    return r / sigma;
}

int Box::makeUnion(int a, int b, int idd[], int size[]) { //union operation in dynamic connectivity
    if((a>=0 && a<npart) && (b>=0 && b<npart)) {
        int i = root(a,idd);
        int j = root(b,idd);
        if(i == j) return 0;
        if(size[i] < size[j]) {
            idd[i] = j;
            size[j] += size[i];
        }
        else {
            idd[j] = i;
            size[i] += size[j];
        }
    }
    else {
        cout << "Error in makeUnion" << endl;
        exit(0);
    }

    return 0;
}

int Box::root(int a,int idd[]) {
    while(a != idd[a]) {
        idd[a] = idd[idd[a]];
        a = idd[a];
    }
    return a;
}

int Box::work(int centralTag, int nn, double* rx, double* ry, double* rz, double* rs, int* cantag, bool& opened) { // function performing voronoi tessellation
//    if(nn < 4) {
//        cout << "less than 4 points in work" << endl;
//        opened = true;
//        return 0;
////        exit(1);
//    }
    
    double ai, bi, ci, di;
    double aj, bj, cj, dj;
    double ak, bk, ck, dk;
    double ab,bc,ca,da,db,dc;
    double alpha[MAXCAN]; 
    double centralSigma, centralSigmasq;
    bool ok = false;
    double det, vxijk, vyijk, vzijk, vrij;
    int iv[MAXVER], jv[MAXVER], kv[MAXVER];
    int edges[MAXCAN];
    
//    ofstream oFile;
//    oFile.open("rij.txt",ios::app);
    
    centralSigma = atom[centralTag].getsigma();
    centralSigmasq = centralSigma * centralSigma;
    
    int v = 0, nv = 0; // number of vertex
    
    for(int i=0; i<=nn-1; i++) { // loop over all candidates
        alpha[i] = 1.0 / 2.0 * (1 + (centralSigmasq - pow(atom[cantag[i]].getsigma(),2.0)) / (4.0 * rs[i]));
        
        atom[cantag[i]].clearVertex();
        atom[cantag[i]].setPointT(rx[i] * alpha[i], ry[i] * alpha[i], rz[i] * alpha[i]);
    }
    
    for(int i=0; i<=nn-3; i++) {
        ai = rx[i];
        bi = ry[i];
        ci = rz[i];
        di = -rs[i]; // negative square distance
        di *= alpha[i];
        for(int j=i+1; j<=nn-2; j++) {
            aj = rx[j];
            bj = ry[j];
            cj = rz[j];
            dj = -rs[j];
            dj *= alpha[j];
            
            ab = ai * bj - aj * bi;
            bc = bi * cj - bj * ci;
            ca = ci * aj - cj * ai;
            da = di * aj - dj * ai;
            db = di * bj - dj * bi;
            dc = di * cj - dj * ci;
            
            for(int k=j+1; k <= nn-1; k++) {
                ak = rx[k];
                bk = ry[k];
                ck = rz[k];
                dk = -rs[k];
                dk *= alpha[k];
                
                det = ak * bc + bk * ca + ck * ab;
                
                if(fabs(det) > TOL) {
                    vxijk = (-dk * bc - ck * db + bk * dc) / det;
                    vyijk = (-dk * ca + ck * da - ak * dc) / det;
                    vzijk = (-dk * ab + ak * db - bk * da) / det;
                    
                    ok = true;
                    int l = 0;
                    while(ok && l <= nn - 1) {
                        if(l != i && l != j && l != k) (rx[l] * vxijk + ry[l] * vyijk + rz[l] * vzijk) < alpha[l] * rs[l] ? ok = true : ok = false;
                        l++;
                    } // if ok == true, include this vertex; if ok == false, exclude it
                    
                    if(ok) {
                        vrij = sqrt(vxijk*vxijk + vyijk*vyijk + vzijk*vzijk);
                        
//                        oFile << vrij / centralSigma << endl;
                        
                        if(vrij < 1.2 * centralSigma) {
                            v++;
                            if(v > MAXVER) {
                                cout << "error: too many vertices" << endl;
                                exit(1);
                            }
                            iv[v - 1] = i;
                            jv[v - 1] = j;
                            kv[v - 1] = k;
                            atom[cantag[i]].push_vertex(vxijk, vyijk, vzijk);
                            atom[cantag[j]].push_vertex(vxijk, vyijk, vzijk);
                            atom[cantag[k]].push_vertex(vxijk, vyijk, vzijk);
                        }
                    }
                }
            }
        }
    }
    
    nv = v;
    
    bool openParticle = false;
    
    if(nv == 0) {
        openParticle = true;
    }
    else {
        for(int i=0;i<=nn-1;i++) {
           if(atom[cantag[i]].isNeighbor()) {
              for(int j=0;j<=nn-1;j++) edges[j] = 0; // all candidates
            
              for(int j=0;j<=nv-1;j++) { // all vertex
                 if(iv[j] == i || jv[j] == i || kv[j] == i) {
                    if(iv[j] != i) edges[iv[j]]++;
                    if(jv[j] != i) edges[jv[j]]++;
                    if(kv[j] != i) edges[kv[j]]++;
                 }
              }
            
              for(int j=0;j<=nn-1;j++)
                  if(edges[j] == 1) {
                     openParticle = true;
                     break;
                  }
            
              if(openParticle) break; // openParticle == true, no enclosing polygon; openParticle == false, being a closed polygon
            }
        }
    }
    
    opened = openParticle;
    
    if(!openParticle) {
        if(nv < 4) {
            cout << "error: less than 4 vertices found in work" << endl;
            cout << "v = " << nv << endl;
            exit(1);
        }
        
        // identify neighboring points
        
        for(int i=0; i<=nn-1; i++) // edges associated with one candidate i
            edges[i] = 0;
        
        for(int i=0; i<=nv-1; i++) {
            edges[iv[i]]++;
            edges[jv[i]]++;
            edges[kv[i]]++;
        }
        
        // check Euler relation
        
        int nf = 0, ne = 0; // number of edges of voronoi polygon
        
        for(int i=0; i<=nn-1; i++)
            if(edges[i] > 0) {
                nf++;
                ne += edges[i];
            }
        
        if(ne % 2 != 0) {
            cout << "Error in work: number of edges is not an even number in Euler relation check" << endl;
            exit(1);
        }
        
        ne /= 2;
        
        if(nv - ne + nf != 2) {
            cout << "Error in work: Euler relation violated" << endl;
            exit(1);
        }
    }
    
//    oFile.close();
    
    return 0;
}

int Box::sort(int centralTag, double* rx, double* ry, double* rz, double* rs, int* nbtag, int& nn) { // sort neighbors into increasing distance order
    bool change = true; // if changed any two items
    int itop = nn - 2; // point to next to top element with array starting from 0, rx(itop) is this element
    int i1;
    double rxi,ryi,rzi,rsi;
    int tagi;
    
    if(nn < 4) return 0;
    
    while(change && itop >= 0) {
        change = false;
        for(int i = 0; i <= itop; i++) {
            i1 = i + 1;
            if(rs[i] > rs[i1]) {
                rxi = rx[i];
                ryi = ry[i];
                rzi = rz[i];
                rsi = rs[i];
                tagi = nbtag[i];
                
                rx[i] = rx[i1];
                ry[i] = ry[i1];
                rz[i] = rz[i1];
                rs[i] = rs[i1];
                nbtag[i] = nbtag[i1];
                
                rx[i1] = rxi;
                ry[i1] = ryi;
                rz[i1] = rzi;
                rs[i1] = rsi;
                nbtag[i1] = tagi;
                
                change = true;
            }
        }
        itop--;  
    }
    
    // the following statements take care of the situation where |sigma1^2 - sigma2^2| < d^2 * 4
    // d -- center to center distance
    int inew = 0, nnnew;
    
    nnnew = nn;
        
    for(int i=0; i<=nn-1; i++) {
        if(fabs(pow(atom[centralTag].getsigma(),2.0) - pow(atom[nbtag[i]].getsigma(),2.0)) >= 4.0 * rs[i]) { // too close, this is invalid candidate
            nnnew--;
            continue;
        }
        else { // normal candidate
            if(inew != i) {
                rx[inew] = rx[i];
                ry[inew] = ry[i];
                rz[inew] = rz[i];
                rs[inew] = rs[i];
                nbtag[inew] = nbtag[i];
            }
            inew++;
        }
    }
    
    nn = nnnew;
    
    return 0;
}
 
double Box::inBoxx(double r) const { // pbc in x direction
    r = r - boxx * round(r / boxx);
    return r;
}

double Box::inBoxy(double r) const { // pbc in y direction
    r = r - boxy * round(r / boxy);
    return r;
}

double Box::inBoxz(double r) const { // pbc in z direction
    r = r - boxz * round(r / boxz);
    return r;
}







