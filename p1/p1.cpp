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
using std::cout;
using std::cerr;
using std::endl;



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
  cout << mol.xyz << endl;

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

  cout << "\nCenter of mass:" << endl;
  auto com = mol.com();
  cout << com.transpose() << endl;
  mol.translate(-com);

  cout << "\nMoment of inertia tensor (amu*bohr^2):" << endl;
  cout << mol.moi_tensor() << endl;
  
  auto pMOI = mol.moi_moms().transpose();
  cout << "\nPrincipal moments of inertia (amu*bohr^2):" << endl;
  cout << pMOI << endl;

  cout << "\nPrincipal moments of inertia (amu*Å^2):" << endl;
  cout << pMOI / angstrom_to_bohr / angstrom_to_bohr << endl;

  cout << "\nPrincipal moments of inertia (g*cm^2):" << endl;
  cout << std::scientific << pMOI * amu_to_g / angstrom_to_bohr / angstrom_to_bohr * 1.0E-16 << endl;

  auto rotconst = mol.rot_const().transpose();
  cout << "\nRotational constants (MHz):" << endl;
  cout << std::scientific << rotconst * speedOfLight * 100 * 1E-06 << endl;

  cout << "\nRotational constants (cm^-1):" << endl;
  cout << std::scientific << rotconst << endl;

/*
  cout << typeid(pMOI).name() << endl;
*/
  return 0;
}

