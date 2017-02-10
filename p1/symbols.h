/* 
 * define global map for element symbol -> atomic number
 * define global vector for atomic number -> element symbol
 * 2017-02-09  vector works properly, map does not
 */
#ifndef _symbols_h_
#define _symbols_h_
#include <vector>
#include <string>
#include <map>
#include <algorithm>
/*
template<typename T, size_t N>
T * end(T (&ra)[N]) {
  return ra + N;
}

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
extern const std::vector<double> atomicMass = {0.0, 1.00794,4.002602,6.941,9.012182,10.811,12.0107,14.0067,15.9994,18.9984032,20.1791,22.98976928,24.3050,26.9815386,28.0855,30.973762,32.065,35.453,39.948,39.0983,40.078,44.955912,47.867,50.9415,51.9961,54.938045,55.845,58.933195,58.6934,63.546,65.38,69.723,72.63,74.92160,78.96,79.904,83.798,85.4678,87.62,88.90585,91.224,92.90638,95.96,98.,101.07,102.90550,106.42,107.8682,112.411,114.818,118.710,121.760,127.60,126.90447,131.293,132.9054519,137.327,138.90547,140.116,140.90765,144.242,145.,150.36,151.964,157.25,158.92535,162.500,164.93032,167.259,168.93421,173.054,174.9668,178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.966569,200.59,204.3833,207.2,208.98040,209.,210.,222.,223.,226.,227.,232.03806,231.03586,238.02891,237.,244.,243.,247.,247.,251.,252.,257.,258.,259.,262.,265.,268.,271.,272.,270.,276.,281.,280.,285.,284.,289.,288.,293.,294.,294.};
extern const std::vector<std::string> elementSymbols = {"X",
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
/* 
const std::map<std::string, int> smap2= {
  {"H", 1},
  {"He", 2},
  {"Li", 3}};
  */
//extern std::map<std::string, int> smap2;
/*
class symbols {
  public:
  static std::map<std::string, int> atomicNumber;
  static std::vector<std::string> elementSymbol;
  
};
std::vector<std::string> symbols::elementSymbol(elementSymbols);
{
for (size_t i = 0; i < elementSymbols.size(); ++i){
  symbols::atomicNumber[elementSymbols[i]] = i;
  symbols::atomicNumber[toUpper(elementSymbols[i])] = i;
  symbols::atomicNumber[toLower(elementSymbols[i])] = i;
}
}
*/
#endif
