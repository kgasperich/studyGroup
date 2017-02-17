#ifndef _molecule_h_
#define _molecule_h_
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
//#include <map>
#include "symbols.h"
#define SW(x) std::setw(x)
using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mxxd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Mx3d;
typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> M33d;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vxd;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vxi;
typedef Eigen::Matrix<double, 3, 1> V3d;
typedef Eigen::SelfAdjointEigenSolver<M33d> EigSol3d;
typedef Eigen::SelfAdjointEigenSolver<Mxxd> EigSolxd;

class Molecule {
  public:
    int natom;
    int znuc;
    int nelec;
    double mass;
    Vxi zvals;
    Vxd masses;
    std::vector<std::string> elnames;
    Mx3d xyz;
    std::string comment;
    Mxxd hess;
    Mxxd hess_mw;
 
    void testfunction();
    V3d rot_const();
    M33d moi_tensor();
    EigSol3d moi_eigs();
    V3d moi_moms();
    int charge();
    void print_xyz();
    void print_zvals();
    void translate(V3d txyz);
    V3d com();
    double bond(int n1, int n2);
    double angle(int n1, int n2, int n3);
    double oop(int n1, int n2, int n3, int n4);
    double dihedral(int n1, int n2, int n3, int n4);
    void read_hessian(const std::string& filename);
    void massweighthessian();
    Molecule(const std::string& filename);
    Molecule();
    ~Molecule();
  private:
    void read_dotxyz(std::istream& is);
};
#endif
