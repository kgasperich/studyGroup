#include "molecule.h"
#include <iostream>
#include <cstdio>

#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
/*
typedef Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> Mx3d;
typedef Eigen::Matrix<double, 3, 1> V3d;


int natom;
int charge;
std::vector<int> zvals;
std::vector<std::string> elnames;
Mx3d xyz;
*/
void Molecule::print_xyz() {
  for (auto i=0; i<natom; i++) {
    printf("%d %3s %8.5f %8.5f %8.5f\n", zvals[i], elnames[i], xyz(i,0), xyz(i,1), xyz(i,2));
  }
}

void Molecule::translate(V3d txyz) {
  for (auto i=0; i<natom; i++) {
    xyz.row(i) += txyz.transpose();
  }
}

Molecule::Molecule(){ }
Molecule::~Molecule(){ }
/*
double bond(int n1, int n2);
double angle(int n1, int n2, int n3);
double oop_angle(int n1, int n2, int n3, int n4);
double dihedral_angle(int n1, int n2, int n3, int n4);

Molecule();
~Molecule();
*/

