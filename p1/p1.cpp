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
// Eigen matrix algebra library
#include <Eigen/Dense>

#define MAXBONDLENGTH 4.0
using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
typedef Eigen::Matrix<double,3,1> V3d;
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

std::vector<Atom> read_geometry(const std::string& filename, bool a2b);

int main(int argc, char *argv[]) {

  bool a2b = false;
  // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
  const auto filename = (argc > 1) ? argv[1] : "acetaldehyde.xyz";
  std::vector<Atom> atoms = read_geometry(filename,a2b);

  // count the number of electrons
  auto nelectron = 0;
  for (auto i = 0; i < atoms.size(); ++i)
    nelectron += atoms[i].atomic_number;
  const auto ndocc = nelectron / 2;

  print_bond_lengths(atoms);
  print_bond_angles(atoms);
  print_oop_angles(atoms);
  print_dihedral_angles(atoms);

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

std::vector<Atom> read_dotxyz(std::istream& is, bool a2b) {
  // line 1 = # of atoms
  size_t natom;
  is >> natom;
  // read off the rest of line 1 and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // line 2 = comment (possibly empty)
  std::string comment;
  std::getline(is, comment);

  std::vector<Atom> atoms(natom);
  for (auto i = 0; i < natom; i++) {
    std::string element_symbol;
    double x, y, z;
    is >> element_symbol >> x >> y >> z;

    // 
    int Z;
    if (isdigit(element_symbol[0])) {
      Z = stoi(element_symbol);
    } else {
    if (element_symbol == "H")
      Z = 1;
    else if (element_symbol == "C")
      Z = 6;
    else if (element_symbol == "O")
      Z = 8;
    else {
      std::cerr << "read_dotxyz: element symbol \"" << element_symbol << "\" is not recognized" << std::endl;
      throw "Did not recognize element symbol in .xyz file";
    }
    }

    atoms[i].atomic_number = Z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr

    atoms[i].xyz << x , y , z;
    const auto angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value
    if(a2b){atoms[i].xyz *= angstrom_to_bohr;}
  }
  return atoms;
}

std::vector<Atom> read_geometry(const std::string& filename, bool a2b) {

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return read_dotxyz(iss, a2b);
  else
    throw "only .xyz files are accepted";
}

void print_bond_lengths(const std::vector<Atom> mol, double rmax) {
  auto enuc = 0.0;
  for (auto i = 0; i < mol.size(); i++)
    for (auto j = 0; j < i; j++) {
      auto rij = mol[i].xyz - mol[j].xyz;
      auto r = rij.norm();
      enuc += mol[i].atomic_number * mol[j].atomic_number / r;
      if (r < rmax) {
        cout << "r" << i << j << " = " << r << endl;
      }
    }
//  cout << "E_nuc = " << enuc << endl;
}

void print_bond_angles(const std::vector<Atom> mol, double rmax) {
  auto N = mol.size();
  for (auto i = 0; i < N; i++)
    for (auto j = 0; j < N; j++) {
      if (j != i) {
// can't use auto here if rji is to be normalized?
//        auto rji = mol[i].xyz - mol[j].xyz;
        V3d rji = mol[i].xyz - mol[j].xyz;
        auto dji = rji.norm();
        if (dji < rmax) {
          rji.normalize();
          for (auto k = 0; k < i; k++) {
            if (k != j) {
//              auto rjk = mol[k].xyz - mol[j].xyz;
              V3d rjk = mol[k].xyz - mol[j].xyz;
              auto djk = rjk.norm();
              if (djk < rmax) {
                rjk.normalize();
                cout << "φ" << i << j << k << " = " << bond_angle(rji,rjk) << endl;
              }
            }
          }
        }
      }
    }
}

void print_oop_angles(const std::vector<Atom> mol, double rmax) {
  auto N = mol.size();
  for (auto k = 0; k < N; k++)
    for (auto i = 0; i < N; i++) {
      if (i != k) {
        V3d rki = mol[i].xyz - mol[k].xyz;
        auto dki = rki.norm();
        if (dki < rmax) {
          rki.normalize();
          for (auto j = 0; j < i; j++) {
            if (j != k) {
              V3d rkj = mol[j].xyz - mol[k].xyz;
              auto dkj = rkj.norm();
              if (dkj < rmax) {
                rkj.normalize();
                for (auto l = 0; l < N; l++) {
                  if ((l != k) && (l != i) && (l != j)) {
                    V3d rkl = mol[l].xyz - mol[k].xyz;
                    auto dkl = rkl.norm();
                    if (dkl < rmax) {
                      rkl.normalize();
                      cout << "θ" << i << j << k << l << " = " << oop_angle(rki,rkj,rkl) << endl;
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


double bond_angle(const Eigen::Vector3d rji, const Eigen::Vector3d rjk, bool normed) {
  auto jidotjk = rji.dot(rjk);
  if (normed) {
    return acos(jidotjk) * 180.0 / M_PI; 
  } else {
    auto jijk = rji.norm()*rjk.norm();
    return acos(jidotjk/jijk) * 180.0 / M_PI; 
  }
}

double oop_angle(const Eigen::Vector3d rki, const Eigen::Vector3d rkj, const Eigen::Vector3d rkl, bool normed) {
  auto kjcrosskl = rkj.cross(rkl);
  kjcrosskl.normalize();
  auto num = kjcrosskl.dot(rki);
  if (normed) {
    return asin(num) * 180.0 / M_PI; 
  } else {
    auto dki = rki.norm();
    return asin(num/dki) * 180.0 / M_PI; 
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
