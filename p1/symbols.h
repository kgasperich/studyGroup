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
