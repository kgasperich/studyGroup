#include "molecule.h"
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <string>
//#include <algorithm>
#include <Eigen/Dense>
#include <cmath>

int Molecule::charge() {
  return znuc - nelec;
}

void Molecule::print_zvals() {
  for (auto i=0; i < natom; i++) {
    cout << SW(3) << zvals[i] << endl;
  }
}
void Molecule::print_xyz() {
  for (auto i=0; i < natom; i++) {
//    printf("%3d %2s %8.5f %8.5f %8.5f\n", zvals[i], elnames[i].c_str(), xyz(i,0), xyz(i,1), xyz(i,2));
    cout 
      << SW(3) << zvals[i]   << " " 
      << SW(2) << elnames[i] << " " 
      << std::fixed << std::setprecision(5)
      << SW(8)  << xyz(i,0)   << " "
      << SW(8)  << xyz(i,1)   << " "
      << SW(8)  << xyz(i,2)
      << endl;
  }
}

M33d Molecule::moi_tensor() {
  M33d I = M33d::Zero();
  auto Irr = (((xyz.rowwise().squaredNorm()).cwiseProduct(masses))).sum();
  for (auto i=0; i < 3; i++) {
    auto Imi = xyz.col(i).cwiseProduct(masses);
    for (auto j=0; j < i; j++) {
      auto Iij = (Imi.cwiseProduct(xyz.col(j))).sum();
      I(i,j) = Iij;
      I(j,i) = Iij;
    }
    auto Iii = (Imi.cwiseProduct(xyz.col(i))).sum();
    I(i,i) = Irr - Iii;
  }
  return I;
}

EigSol3d Molecule::moi_eigs() {
  M33d I = moi_tensor();
  return EigSol3d(I);
}

V3d Molecule::moi_moms() {
  M33d I = moi_tensor();
  return EigSol3d(I).eigenvalues();
}

//in cm^-1
V3d Molecule::rot_const() {
  V3d Ii = moi_moms();
  return Ii.array().inverse() * planckConstantJs / (8.0 * M_PI * M_PI * speedOfLight) /(amu_to_kg ) * (m_to_bohr * m_to_bohr) / 100.0;
}

void Molecule::translate(V3d txyz) {
  for (auto i=0; i < natom; i++) {
    xyz.row(i) += txyz.transpose();
  }
}

double Molecule::bond(int n1, int n2) {
  return (xyz.row(n2)-xyz.row(n1)).norm();
}

// angle 1-2-3
double Molecule::angle(int n1, int n2, int n3) {
  auto r21 = xyz.row(n1) - xyz.row(n2);
  auto r23 = xyz.row(n3) - xyz.row(n2);
  auto costheta = r21.dot(r23)/(r21.norm()*r23.norm());
  return acos(costheta) * 180.0 / M_PI; 
}

// angle between r31 and plane 2-3-4
double Molecule::oop(int n1, int n2, int n3, int n4) {
  auto r32 = (xyz.row(n2) - xyz.row(n3)).transpose();
  auto r34 = (xyz.row(n4) - xyz.row(n3)).transpose();
  V3d r31 = (xyz.row(n1) - xyz.row(n3)).transpose();
  auto s234 = r32.cross(r34);
  s234.normalize();
  r31.normalize();
  return asin(s234.dot(r31)) * 180.0 / M_PI;
}

double Molecule::dihedral(int n1, int n2, int n3, int n4) {
  auto r21 = (xyz.row(n1) - xyz.row(n2)).transpose();
  auto r23 = (xyz.row(n3) - xyz.row(n2)).transpose();
  auto r34 = (xyz.row(n4) - xyz.row(n3)).transpose();
  auto s123 = r21.cross(r23);
  auto s234 = r34.cross(r23);
  s123.normalize();
  s234.normalize();
  auto sign = 1;
  if ((s123.cross(s234)).dot(r23) < 0.0) sign = -1;
  return sign * acos(s123.dot(s234)) * 180.0 / M_PI;
}

V3d Molecule::com() {
  V3d xyzcom = V3d::Zero();
  for (auto i=0; i < natom; i++) {
    xyzcom += xyz.row(i) * masses[i];
  }
  return xyzcom/mass;
}
//Molecule::Molecule(int n, int q): 
//  natom(n), charge(q), 
//  zvals(n,0), masses(n,0), elnames(n), xyz(n,3) { }

Molecule::Molecule(){ }

Molecule::~Molecule(){ }

void Molecule::read_dotxyz(std::istream& is) {
  // line 1 = # of atoms
  int l1;
  is >> l1;
  bool a2b;
  if (l1 < 0) {
    a2b = false;
    natom = -l1;
    cout << "Input is in atomic units" << endl;
  } else {
    a2b = true;
    natom = l1;
    cout << "Input is in angstroms" << endl;
  }
  zvals.resize(natom);
  elnames.resize(natom);
  masses.resize(natom);
  xyz.resize(natom,3);
  
  // read off the rest of line 1 and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // line 2 = comment (possibly empty)
  std::getline(is, comment);

  for (auto i = 0; i < natom; i++) {
    std::string element_symbol;
    double x, y, z;
    is >> element_symbol >> x >> y >> z;

    // if given as numeric value, convert to int
    if (isdigit(element_symbol[0])) {
      zvals[i] = stoi(element_symbol);
    } else {
    // if given as symbol, use map to get value
      try {
        zvals[i] = symbolToZ.at(element_symbol);
      } catch (std::exception& oe) {
        cerr << "\nError: element " << element_symbol << " not recognized\n" << oe.what() << endl;
        throw "Did not recognize element symbol in .xyz file";
      }
    }
    elnames[i] = elementSymbols[zvals[i]];
    masses[i] = atomicMassIso[zvals[i]];
    mass += masses[i];
    znuc += zvals[i];
    nelec += zvals[i];
    xyz.row(i) << x , y , z;
  }
  if(a2b){xyz *= angstrom_to_bohr;}
}

Molecule::Molecule(const std::string& filename): mass(0), znuc(0), nelec(0) {

  std::cout << "Reading geometry from " << filename << std::endl;
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
    read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

