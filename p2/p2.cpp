/* 
 * Crawford exercise 2
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <typeinfo>

#include "symbols.h"
#include "molecule.h"

#include <Eigen/Dense>

#define RMAX 4.0
using std::cout;
using std::cerr;
using std::endl;



int main(int argc, char *argv[]) {

  const auto filename = (argc > 1) ? argv[1] : "h2o_geom.xyz";
  const auto hessfilename = (argc > 2) ? argv[2] : "h2o_hessian.txt";

  Molecule mol(filename);
  mol.read_hessian(hessfilename);

  cout << "\ngeometry:" << endl;
  mol.print_xyz();

  cout << "\nhessian:" << endl;
  cout << mol.hess << endl;

  cout << "\nmass-weighted hessian:" << endl;
  cout << mol.hess_mw << endl;

  EigSolxd hessEigSol(mol.hess_mw);
  auto hessEigVals = hessEigSol.eigenvalues();
  cout << "\neigvals (Eh/(bohr^2 * amu)):" << endl;
  cout << hessEigVals << endl;

  auto fConv = sqrt(hartree_to_J * m_to_bohr * m_to_bohr / amu_to_kg)/(2.0 * M_PI * speedOfLight * 100);
  auto freq = hessEigVals.cwiseSqrt() * fConv;
  cout << "\nfrequencies (cm^-1):" << endl;
  cout << freq << endl;

  return 0;
}

