#ifndef DTFRPARSER_H_
#define DTFRPARSER_H_

#include <vector>
#include <string>

#include "nuclidedata.h"

class DtfrParser
{
public:
    DtfrParser();
    ~DtfrParser();

    bool openFile(std::string filename);
    bool parseFile(std::string filename);
    bool parseHeader();
    bool parseData();
    void closeFile();

    inline const std::string& getFilename() const { return m_filename; }
    inline std::vector<int> getZaids() const { return zAids; }
    inline int getNumberNuclides() const { return nNuclides; }
    inline int getPnOrder() const { return nPn; }
    inline int getGammaEnergyGroups() const { return gGroups; }
    inline const std::vector<float>& getGammaEnergy() const { return gBounds; }
    inline const NuclideData* getNuclideEntry(const int indx) const
                              {
                                  return data[indx*(nPn+1)];
                              }

    int getIndexbyZaid(int zaid) const;
    void debugGammaEnergies();
    NuclideData* getData(unsigned int indx, unsigned int ln) const;

protected:
    std::ifstream binfile;           //cross section input file
    std::string m_filename;          //cross section data file name
    int gGroups;                     //number of energy groups
    int nNuclides;                   //number of nuclides in DTFR data file
    int nPn;                         //order of Legendre expansion start from zero
    std::vector<float> gBounds;      //energy bin boundary
    std::vector<int> zAids;          //atomic number
    std::vector<NuclideData*> data;  //cross section data including 1D,2D scatter matrix

};
#endif // DTFRPARSER_H_
