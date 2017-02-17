/* 
 * Crawford exercise 2
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
  return 0;
}

