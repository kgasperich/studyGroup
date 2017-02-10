#ifndef _molecule_h_
#define _molecule_h_
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Mx3d;
typedef Eigen::Matrix<double, 3, 1> V3d;

class Molecule {
  public:
    int natom;
    int charge;
    std::vector<int> zvals;
    std::vector<std::string> elnames;
    Mx3d xyz;

    void print_xyz();
    void translate(V3d txyz);
    double bond(int n1, int n2);
    double angle(int n1, int n2, int n3);
    double oop_angle(int n1, int n2, int n3, int n4);
    double dihedral_angle(int n1, int n2, int n3, int n4);

    Molecule();
    ~Molecule();
};
#endif
