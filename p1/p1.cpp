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

// Eigen matrix algebra library
#include <Eigen/Dense>

using std::cout;
using std::cerr;
using std::endl;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library

struct Atom {
    int atomic_number;
    Eigen::Vector3d xyz;
};

void print_bond_lengths(const std::vector<Atom> mol);
void print_bond_angles(const std::vector<Atom> mol, double rmax = 3.0);
double angle(const Eigen::Vector3d ji, const Eigen::Vector3d jk);
std::vector<Atom> read_geometry(const std::string& filename, bool a2b);
int main(int argc, char *argv[]) {

  bool a2b = false;
    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename,a2b);
    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;
    print_bond_lengths(atoms);
    print_bond_angles(atoms);


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

void print_bond_lengths(const std::vector<Atom> mol) {
  auto enuc = 0.0;
  for (auto i = 0; i < mol.size(); i++)
    for (auto j = 0; j < i; j++) {
      auto rij = mol[i].xyz - mol[j].xyz;
      auto r = rij.norm();
      enuc += mol[i].atomic_number * mol[j].atomic_number / r;
      cout << "r" << i << j << " = " << r << endl;
    }
//  cout << "E_nuc = " << enuc << endl;
}

void print_bond_angles(const std::vector<Atom> mol, double rmax) {
  for (auto i = 0; i < mol.size(); i++)
    for (auto j = 0; j < i; j++) {
      auto rij = mol[j].xyz - mol[i].xyz;
      auto dij = rij.norm();
      for (auto k = 0; k < j; k++) {
        auto rjk = mol[k].xyz - mol[j].xyz;
        auto rki = mol[i].xyz - mol[k].xyz;
        auto djk = rjk.norm();
        auto dki = rki.norm();
        if (dij < rmax) {
          if (djk < rmax) {
            cout << "θ" << i << j << k << " = " << angle(-1.0*rij,rjk) << endl;
          }
          if (dki < rmax) {
            cout << "θ" << j << i << k << " = " << angle(rij,-1.0*rki) << endl;
          }
        } else if ((djk < rmax) && (dki < rmax)) {
            cout << "θ" << i << k << j << " = " << angle(rki,-1.0*rjk) << endl;
        }

      }
    }
}

double angle(const Eigen::Vector3d ji, const Eigen::Vector3d jk) {
  auto jidotjk = ji.dot(jk);
  auto jijk = ji.norm()*jk.norm();
  auto tijk = acos(jidotjk/jijk) * 180.0 / M_PI; 
  return tijk;
}
