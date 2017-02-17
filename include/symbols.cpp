/* 
 * define global map for element symbol -> atomic number
 * define global vector for atomic number -> element symbol
 * 2017-02-09  vector works properly, map does not
 */
#include "symbols.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
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
// mass of most common isotope (from Wolfram IsotopeData[])
const std::vector<double> atomicMassIso = {0.000,1.00782503207,4.00260325415,7.016004548,9.012182201,11.009305406,12.,14.00307400478,15.99491461956,18.998403224,19.99244017542,22.98976928087,23.985041699,26.981538627,27.97692653246,30.973761629,31.972070999,34.968852682,39.96238312251,38.963706679,39.962590983,44.955911909,47.947946281,50.943959507,51.940507472,54.938045141,55.934937475,58.933195048,57.935342907,62.929597474,63.929142222,68.925573587,73.921177767,74.921596478,79.916521271,78.918337087,83.911506687,84.911789737,87.905612124,88.905848295,89.904704416,92.906378058,97.905408169,117.951480,101.904349312,102.905504292,105.903485715,106.90509682,113.90335854,114.903878484,119.902194676,120.903815686,129.906224399,126.904472681,131.904153457,132.905451932,137.905247237,138.906353267,139.905438706,140.907652769,141.907723297,162.953680,151.919732425,152.921230339,157.924103912,158.925346757,163.929174751,164.93032207,165.930293061,168.93421325,173.938862089,174.940771819,179.946549953,180.947995763,183.950931188,186.955753109,191.96148069,192.96292643,194.964791134,196.966568662,201.970643011,204.974427541,207.976652071,208.980398734,220.016602,223.025190,228.037986,232.049772,234.050704,236.055296,232.038055325,240.060980,238.050788247,244.067850,247.074070,249.078480,252.084870,254.090600,256.093440,258.099520,260.102678,262.108865,264.112345,266.119305,268.123644,270.130712,273.138220,275.144250,277.149841,279.156193,281.162061,283.168415,285.174105,287.181045,289.187279,291.194384,292.199786,292.207549,293.214670};
// average mass (from Wolfram ElementData[])
const std::vector<double> atomicMassAvg = {0.0,1.00794,4.002602,6.941,9.012182,10.811,12.0107,14.0067,15.9994,18.9984032,20.1791,22.98976928,24.3050,26.9815386,28.0855,30.973762,32.065,35.453,39.948,39.0983,40.078,44.955912,47.867,50.9415,51.9961,54.938045,55.845,58.933195,58.6934,63.546,65.38,69.723,72.63,74.92160,78.96,79.904,83.798,85.4678,87.62,88.90585,91.224,92.90638,95.96,98.,101.07,102.90550,106.42,107.8682,112.411,114.818,118.710,121.760,127.60,126.90447,131.293,132.9054519,137.327,138.90547,140.116,140.90765,144.242,145.,150.36,151.964,157.25,158.92535,162.500,164.93032,167.259,168.93421,173.054,174.9668,178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.966569,200.59,204.3833,207.2,208.98040,209.,210.,222.,223.,226.,227.,232.03806,231.03586,238.02891,237.,244.,243.,247.,247.,251.,252.,257.,258.,259.,262.,265.,268.,271.,272.,270.,276.,281.,280.,285.,284.,289.,288.,293.,294.,294.};
const std::vector<std::string> elementSymbols = {"X",
//extern const char *elementSymbols[] = {"X",
    "H",                                                                                 "He",
    "Li","Be",                                                  "B", "C", "N", "O", "F", "Ne",
    "Na","Mg",                                                  "Al","Si","P", "S", "Cl","Ar",
    "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I", "Xe",
    "Cs","Ba",
              "La","Ce","Pr","Nd","Pm","Sm","Eu","Gc","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                   "Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra",
              "Ac","Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
                   "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};


int symbolToZ(std::string symbol) {
  auto i = std::distance(elementSymbols.begin(), 
      std::find(elementSymbols.begin(), elementSymbols.end(), symbol));
  if ( i == elementSymbols.size()) {
    std::cerr << "Element symbol \"" << symbol << "\" not recognized" << std::endl;
    throw std::runtime_error("error reading element symbol.");
  }
  return i;
}