#ifndef NUCLIDEDATA_H_
#define NUCLIDEDATA_H_

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

class NuclideData
{
public:
    NuclideData();
    ~NuclideData();

    bool parse(std::ifstream &binfile, int gGroups);

    inline const std::vector<float> &getGammaScatter1D() const {return gScatter1d;}
    inline const std::vector<float> &getGammaScatter2D() const {return gScatter2d;}
    inline const std::vector<float> &getGammaXs() const {return tolxsection; }
    inline const std::vector<float> &getGammaCohScatter() const {return cohxsection;}
    inline const std::vector<float> &getGammaIncScatter() const {return incxsection;}
    inline const std::vector<float> &getGammaPhotoelec() const {return phexsection;}
    inline const std::vector<float> &getGammaPairprod() const {return parxsection; }

protected:
    std::vector<float> gScatter1d;    //group scatter 1D matrix
    std::vector<float> gScatter2d;    //group scatter 2D matrix
    std::vector<float> tolxsection;   //total cross section
    std::vector<float> cohxsection;   //coherent scatter cross section
    std::vector<float> incxsection;   //incoherent scatter cross section
    std::vector<float> phexsection;   //photo-electric cross section
    std::vector<float> parxsection;   //pair-production cross section
};


#endif // NUCLIDEDATA_H_
