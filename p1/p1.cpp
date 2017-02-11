/* Crawford exercise 1
 *
 * reused code from github.com/evaleev/libint/tests/hartree-fock/hartree-fock.cc 
 *   for reading in .xyz files, counting electrons, and calculating E_nuc
 *
 *
 * 2017-02-09 framework in place to read .xyz, calculate bond lengths, calculate E_nuc
 * TODO: handle angstrom/bohr conversion properly
 * TODO: make global map for mapping from element symbols to numbers {{"H",1},{"He",2},...}
 * TODO: make global vector for mapping from numbers to symbols {"X","H","He","Li",...}
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <typeinfo>
//#include "symbols.h"
#include "molecule.h"
// Eigen matrix algebra library
#include <Eigen/Dense>
#define RMAX 4.0
#define MAXBONDLENGTH 3.0
using std::cout;
using std::cerr;
using std::endl;

//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
//typedef Eigen::Matrix<double,3,1> V3d;
struct Atom {
    int atomic_number;
    Eigen::Vector3d xyz;
};

void print_bond_lengths(const std::vector<Atom> mol, double rmax = MAXBONDLENGTH);
void print_bond_angles(const std::vector<Atom> mol, double rmax = MAXBONDLENGTH);
void print_oop_angles(const std::vector<Atom> mol, double rmax = MAXBONDLENGTH);
void print_dihedral_angles(const std::vector<Atom> mol, double rmax = MAXBONDLENGTH);
double bond_angle(const Eigen::Vector3d rji, const Eigen::Vector3d rjk, bool normed = true);
double oop_angle(const Eigen::Vector3d rki, const Eigen::Vector3d rkj, const Eigen::Vector3d rkl, bool normed = true);
double dihedral_angle(const Eigen::Vector3d rji,const Eigen::Vector3d rjk,const Eigen::Vector3d rkl);
V3d center_of_mass(const std::vector<Atom> mol);
void translate(std::vector<Atom> &mol, V3d txyz);
void print_xyz(const std::vector<Atom> mol);
std::vector<Atom> read_geometry(const std::string& filename, bool a2b);

int main(int argc, char *argv[]) {

  // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
  const auto filename = (argc > 1) ? argv[1] : "acetaldehyde.xyz";

  // count the number of electrons
//  auto nelectron = 0;
//  for (auto i = 0; i < atoms.size(); ++i)
//    nelectron += atoms[i].atomic_number;
//  const auto ndocc = nelectron / 2;
  Molecule mol(filename);
  mol.print_xyz();

  cout << "\nBond lengths:" << endl;
  for (auto i=0; i < mol.natom; i++) {
    for (auto j=0; j < i; j++) {
      cout << "r " << i << " " << j << "  " << mol.bond(i,j) << endl;
    }
  }

  cout << "\nBond angles (correct):" << endl;
  for (auto i = 0; i < mol.natom; i++)
    for (auto j = 0; j < mol.natom; j++) {
      if (j != i) {
        if (mol.bond(i,j) < RMAX) {
          for (auto k = 0; k < i; k++) {
            if (k != j) {
              if (mol.bond(j,k) < RMAX) {
                cout << "φ " << i << " " << j << " "<< k << "  " << mol.angle(i,j,k) << endl;
              }
            }
          }
        }
      }
    }
  cout << "\nBond angles (Crawford):" << endl;
  for (auto i = 0; i < mol.natom; i++)
    for (auto j = 0; j < i; j++) {
      if (mol.bond(i,j) < 4.0) {
        for (auto k = 0; k < j; k++) {
          if (mol.bond(j,k) < 4.0) {
            cout << "φ " << i << " " << j << " "<< k << "  " << mol.angle(i,j,k) << endl;
          }
        }
      }
    }

  cout << "\nOOP angles:" << endl;
  for (auto k = 0; k < mol.natom; k++)
    for (auto i = 0; i < mol.natom; i++) {
      if (i != k) {
        if (mol.bond(i,k) < RMAX) {
          for (auto j = 0; j < i; j++) {
            if (j != k) {
              if (mol.bond(j,k) < RMAX) {
                for (auto l = 0; l < mol.natom; l++) {
                  if ((l != k) && (l != i) && (l != j)) {
                    if (mol.bond(k,l) < RMAX) {
                      cout << "θ " 
                        << i << " " << j << " " << k << " " << l << "  " 
                        << mol.oop(i,j,k,l) << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  cout << "\nOOP angles (Crawford):" << endl;
  for (auto i = 0; i < mol.natom; i++)
    for (auto k = 0; k < mol.natom; k++) {
      if (i != k) {
        if (mol.bond(i,k) < RMAX) {
          for (auto j = 0; j < mol.natom; j++) {
            if ((j != k) && (j != i)) {
              if (mol.bond(j,k) < RMAX) {
                for (auto l = 0; l < j; l++) {
                  if ((l != k) && (l != i)) {
                    if (mol.bond(k,l) < RMAX) {
                      cout << "θ " 
                        << i << " " << j << " " << k << " " << l << "  " 
                        << SW(10) << mol.oop(i,j,k,l) << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  cout << "\ndihedral angles:" << endl;
  for (auto j = 0; j < mol.natom; j++)
    for (auto k = 0; k < j; k++) {
      if (k != j) {
        if (mol.bond(j,k) < RMAX) {
          for (auto i = 0; i < mol.natom; i++) {
            if ((i != j) && (i != k)) {
              if (mol.bond(j,i) < RMAX) {
                for (auto l = 0; l < mol.natom; l++) {
                  if ((l != k) && (l != i) && (l != j)) {
                    if (mol.bond(k,l) < RMAX) {
                      cout << "τ "
                        << i << " " 
                        << j << " " 
                        << k << " " 
                        << l << "  " << SW(10) << mol.dihedral(i,j,k,l) << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  cout << "\ndihedral angles (Crawford):" << endl;
  for (auto i = 0; i < mol.natom; i++)
    for (auto j = 0; j < i; j++) {
      if (mol.bond(j,i) < RMAX) {
        for (auto k = 0; k < j; k++) {
          if (mol.bond(j,k) < RMAX) {
            for (auto l = 0; l < k; l++) {
              if (mol.bond(k,l) < RMAX) {
                cout << "τ "
                  << i << " " 
                  << j << " " 
                  << k << " " 
                  << l << "  " << SW(10) << mol.dihedral(i,j,k,l) << endl;
              }
            }
          }
        }
      }
    }

  cout << "\ncenter of mass:" << endl;
  auto com = mol.com();
  cout << com.transpose() << endl;
  mol.translate(-com);
  com = mol.com();
  cout << com.transpose() << endl;

//  auto r1 = mol.xyz.row(1);
//  cout << r1 << endl;
//  cout << typeid(r1).name() << endl;
//  cout << r1.norm() << endl;
/*  
  cout << atoms[5].xyz.array() * atoms[6].xyz.array() << endl;
  print_bond_lengths(atoms);
  print_bond_angles(atoms);
  print_oop_angles(atoms);
  print_dihedral_angles(atoms);
  print_xyz(atoms);
  auto txyz = center_of_mass(atoms);
  cout << txyz.transpose() << endl;
  translate(atoms, -1*txyz);
  print_xyz(atoms);
*/
/*
  auto r1 = atoms[1].xyz;
  V3d r1b = atoms[1].xyz;
  Eigen::Vector3d r1c = atoms[1].xyz;

  cout << typeid(r1).name() << endl;
  cout << typeid(r1b).name() << endl;
  cout << typeid(r1c).name() << endl;
*/
  return 0;
}


void print_dihedral_angles(const std::vector<Atom> mol, double rmax) {
  auto N = mol.size();
  for (auto j = 0; j < N; j++)
    for (auto k = 0; k < j; k++) {
      if (k != j) {
        V3d rjk = mol[k].xyz - mol[j].xyz;
        auto djk = rjk.norm();
        if (djk < rmax) {
          for (auto i = 0; i < N; i++) {
            if ((i != j) && (i != k)) {
              V3d rji = mol[i].xyz - mol[j].xyz;
              auto dji = rji.norm();
              if (dji < rmax) {
                for (auto l = 0; l < N; l++) {
                  if ((l != k) && (l != i) && (l != j)) {
                    V3d rkl = mol[l].xyz - mol[k].xyz;
                    auto dkl = rkl.norm();
                    if (dkl < rmax) {
                      cout << "τ" << i << j << k << l << " = " << dihedral_angle(rji,rjk,rkl) << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
}

V3d center_of_mass(const std::vector<Atom> mol) {
  V3d xyzcom = V3d::Zero();
  auto mtot = 0.0;
  for (auto i = 0; i < mol.size(); i++) {
    auto mi = atomicMassIso[mol[i].atomic_number];
    xyzcom += mi * mol[i].xyz;
    mtot += mi;
  }
  return xyzcom / mtot;
}

void translate(std::vector<Atom> &mol, V3d txyz) {
  for (auto i = 0; i < mol.size(); i++) {
    mol[i].xyz += txyz;
  }
}

double dihedral_angle(const Eigen::Vector3d rji,const Eigen::Vector3d rjk,const Eigen::Vector3d rkl) {
  //TODO: compare n1 x n2 to rjk to determine sign of angle
  V3d n1 = rji.cross(rjk);
  n1.normalize();
  V3d n2 = -1.0 * rjk.cross(rkl);
  n2.normalize();
  return acos(n1.dot(n2)) * 180.0 / M_PI;
}
