/* 
 * define global map for element symbol -> atomic number
 * define global vector for atomic number -> element symbol
 * 2017-02-09  vector works properly, map does not
 */
#ifndef _symbols_h_
#define _symbols_h_
#include <vector>
#include <string>
/*
std::string toUpper(std::string s) {
  std::string s2(s);
  std::transform(s2.begin(), s2.end(), s2.begin(), ::toupper);
  return s2;
}

std::string toLower(std::string s) {
  std::string s2(s);
  std::transform(s2.begin(), s2.end(), s2.begin(), ::tolower);
  return s2;
}
*/
const auto planckConstanteVs = 4.135667662E-15; // 2014 CODATA value (eV s)
const auto planckConstantJs = 6.62607004E-34; // 2014 CODATA value (J s)
const auto speedOfLight = 2.99792458E+08; // 2014 CODATA value (m/s)
const auto angstrom_to_bohr = 1 / 0.52917721067; // 2014 CODATA value
const auto m_to_bohr = angstrom_to_bohr * 1.0E+10; // 2014 CODATA value
const auto amu_to_g = 1.66053904E-24; // 2014 CODATA value
const auto amu_to_kg = amu_to_g * 1.0E-03; // 2014 CODATA value
const auto Me_to_kg = 9.10938356E-31; // 2014 CODATA value
const auto hartree_to_J = 4.35974465E-18; // 2014 CODATA value
// mass of most common isotope (from Wolfram IsotopeData[])
extern const std::vector<double> atomicMassIso;
// average mass (from Wolfram ElementData[])
extern const std::vector<double> atomicMassAvg;
extern const std::vector<std::string> elementSymbols;
int symbolToZ(std::string symbol);

#endif
